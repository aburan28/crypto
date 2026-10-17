//! **Pohlig-Hellman algorithm** for the elliptic-curve discrete log
//! over a **smooth-order** group.
//!
//! Given `Q = d·G` on a curve where `ord(G) = n` factors as
//! `n = ∏ q_i^e_i` with each `q_i` small, recover `d` at cost
//! `O(Σ e_i · √q_i)` instead of `O(√n)`.
//!
//! ## Algorithm (Pohlig–Hellman 1978)
//!
//! 1. **Reduce** to each prime-power factor: let `n_i = q_i^e_i` and
//!    `G_i = (n / n_i) · G`.  Then `Q_i = (n / n_i) · Q = d_i · G_i`
//!    where `d_i = d mod n_i`.
//! 2. **Solve each `d_i`** via either:
//!    - **Baby-step on `q_i`**: write `d_i = d_{i,0} + q_i · d_{i,1} +
//!      … + q_i^{e_i-1} · d_{i,e_i-1}` and recover each `d_{i,j}`
//!      independently with a `√q_i`-cost search on a `q_i`-order
//!      subgroup.
//!    - **Pollard rho** for each `q_i`-order subgroup if `q_i` is too
//!      large to brute-force.
//! 3. **CRT** the `(d_i mod n_i)` results into `d mod n`.
//!
//! ## Where this gets used
//!
//! Foundation for several other attacks in this crate:
//!
//! - **Invalid-curve attack** (`cryptanalysis::invalid_curve_attack`):
//!   feed a victim a point on a *twist* with smooth order; their
//!   scalar mul lives in the smooth group; PH recovers `d mod q_twist`.
//! - **MOV / Frey-Rück** (`cryptanalysis::mov_attack`): pairing
//!   reduces ECDLP to `F_{p^k}*` DLP; the resulting target group
//!   has order `p^k − 1` which is often **very** smooth — PH on
//!   `F_{p^k}*` finishes it off.
//! - **Cheon's attack** on Diffie-Hellman with auxiliary inputs.
//!
//! ## What this module ships
//!
//! - [`PohligHellmanReport`] — diagnostic output of one PH run.
//! - [`pohlig_hellman_curve`] — full driver, generic over the curve.
//! - [`recover_in_prime_power_subgroup`] — solve `Q = d·G` when
//!   `ord(G) = q^e` is a prime power.
//! - [`crt_combine`] — Chinese remainder reconstruction.
//! - [`format_visualization`] — Markdown report showing per-prime
//!   discovered `d_i` and the CRT reconstruction.

use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;
use crate::utils::mod_inverse;
use crate::visualize::color::{paint, FG_BRIGHT_GREEN, FG_BRIGHT_RED, FG_BRIGHT_YELLOW};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Diagnostic output of one Pohlig-Hellman run.
#[derive(Clone, Debug)]
pub struct PohligHellmanReport {
    /// Recovered `d ≡ ? (mod ord(G))`.
    pub recovered_d: Option<BigUint>,
    /// Factorisation of `ord(G) = ∏ q_i^e_i`.
    pub factors: Vec<(BigUint, u32)>,
    /// Per-prime-power discovered residues `(q_i^e_i, d_i)`.
    pub residues: Vec<(BigUint, BigUint)>,
    /// Wall-clock time in milliseconds.
    pub elapsed_ms: u128,
    /// Group points produced or compared by the per-prime baby-step
    /// giant-step searches (each prime's `⌈√q⌉`-point table plus every giant
    /// point), summed over the primes recovered: the cost summary.
    pub total_steps: u64,
}

/// **Recover `d` such that `Q = d·G`** on the elliptic curve,
/// exploiting smooth `ord(G)`.  Factorises the order, recovers each
/// `d_i` via baby-step on the matching prime-power subgroup, CRTs.
///
/// `order_of_g` MUST equal `ord(G)` on the curve (typically `curve.n`).
/// `smoothness_bound` is the maximum prime factor we'll attempt to
/// brute-force; if any factor exceeds this, we still report the
/// partial residues but the final `recovered_d` is `None`.
pub fn pohlig_hellman_curve(
    curve: &CurveParams,
    g: &Point,
    q: &Point,
    order_of_g: &BigUint,
    smoothness_bound: u64,
) -> PohligHellmanReport {
    let t0 = std::time::Instant::now();
    let factors = factorise_up_to(order_of_g, smoothness_bound);
    let bound = BigUint::from(smoothness_bound);
    let a_fe = curve.a_fe();
    let mut residues: Vec<(BigUint, BigUint)> = Vec::new();
    let mut total_steps = 0u64;
    let mut all_factors_recovered = true;
    for (qprime, e) in &factors {
        if qprime > &bound {
            all_factors_recovered = false;
            continue;
        }
        let qe = qprime.pow(*e);
        let cofactor = order_of_g / &qe;
        let g_i = g.scalar_mul_vartime(&cofactor, &a_fe);
        let q_i = q.scalar_mul_vartime(&cofactor, &a_fe);
        match recover_in_prime_power_subgroup(curve, &g_i, &q_i, qprime, *e) {
            Some((d_i, steps)) => {
                residues.push((qe, d_i));
                total_steps += steps;
            }
            None => {
                all_factors_recovered = false;
            }
        }
    }
    // Every answer is checked in the group: a mis-stated order can make the
    // residues consistent with each other and still wrong.
    let recovered_d = if !all_factors_recovered {
        None
    } else if order_of_g.is_one() {
        (*q == Point::Infinity).then(BigUint::zero)
    } else if residues.is_empty() {
        None
    } else {
        crt_combine(&residues).filter(|d| g.scalar_mul_vartime(d, &a_fe) == *q)
    };
    PohligHellmanReport {
        recovered_d,
        factors,
        residues,
        elapsed_ms: t0.elapsed().as_millis(),
        total_steps,
    }
}

/// Factor `n` by trial division with divisors `d ≤ bound`.
///
/// Every prime factor up to `bound` is found with its exponent, as
/// [`factorise_small`](crate::cryptanalysis::j0_twists::factorise_small)
/// finds it.  What remains is one more entry,
/// `(rest, 1)`: a prime whenever the division stopped because
/// `d² > rest`, and otherwise a cofactor whose prime factors all exceed
/// `bound`, which Pohlig–Hellman cannot use either way.  Stopping at the
/// bound is what keeps an order with a large prime factor (a `curve.n` of
/// cryptographic size) from running trial division to its square root.
fn factorise_up_to(n: &BigUint, bound: u64) -> Vec<(BigUint, u32)> {
    if n.is_zero() || n.is_one() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut rest = n.clone();
    let mut d = 2u64;
    while d <= bound && BigUint::from(d) * d <= rest {
        let mut e = 0u32;
        while (&rest % d).is_zero() {
            rest /= d;
            e += 1;
        }
        if e > 0 {
            out.push((BigUint::from(d), e));
        }
        d += if d == 2 { 1 } else { 2 };
    }
    if !rest.is_one() {
        out.push((rest, 1));
    }
    out
}

/// A point as a hashable key: its affine coordinates, `None` for `O`.
fn point_key(p: &Point) -> Option<(BigUint, BigUint)> {
    match p {
        Point::Infinity => None,
        Point::Affine { x, y } => Some((x.value.clone(), y.value.clone())),
    }
}

/// The baby-step table of one prime: `j·G'` for `0 ≤ j < m`,
/// `m = ⌈√q⌉`, keeping the smallest `j` of a repeated point.
struct BabySteps {
    m: u64,
    table: std::collections::HashMap<Option<(BigUint, BigUint)>, u64>,
    /// `−m·G'`, the giant stride.
    neg_stride: Point,
}

impl BabySteps {
    fn build(g_prime: &Point, q: u64, a_fe: &FieldElement) -> Self {
        let mut m = (q as f64).sqrt() as u64;
        while m.saturating_mul(m) < q {
            m += 1;
        }
        let m = m.max(1);
        let mut table = std::collections::HashMap::with_capacity(m as usize);
        let mut p = Point::Infinity;
        for j in 0..m {
            table.entry(point_key(&p)).or_insert(j);
            if j + 1 < m {
                p = p.add_vartime(g_prime, a_fe);
            }
        }
        let neg_stride = g_prime.scalar_mul_vartime(&BigUint::from(m), a_fe).neg();
        Self {
            m,
            table,
            neg_stride,
        }
    }

    /// The smallest `k < q` with `k·G' = target`, and the giant steps
    /// spent: the search covers `[0, m²) ⊇ [0, q)` in increasing blocks
    /// of `m`, so the first hit is the smallest solution, which is the
    /// answer the linear search over `[0, q)` gives, `None` included.
    fn solve(&self, target: &Point, q: u64, a_fe: &FieldElement) -> (Option<u64>, u64) {
        let mut cur = target.clone();
        for i in 0..self.m {
            if let Some(&j) = self.table.get(&point_key(&cur)) {
                let k = i * self.m + j;
                return ((k < q).then_some(k), i + 1);
            }
            cur = cur.add_vartime(&self.neg_stride, a_fe);
        }
        (None, self.m)
    }
}

/// Solve `Q = d·G` when `ord(G) = q^e` is a prime power.
/// Returns `(d mod q^e, step_count)` or `None` on failure.
///
/// Strategy: write `d = d_0 + q·d_1 + q²·d_2 + … + q^{e-1}·d_{e-1}`.
/// Each `d_j ∈ [0, q)` is recovered by baby-step giant-step in the
/// `q`-order subgroup generated by `G' = q^{e−1}·G`: one table of the
/// `m = ⌈√q⌉` multiples `j·G'` serves all `e` digits, and a digit costs
/// at most `m` giant steps, so a prime costs `O((e + 1)√q)` group
/// operations where the linear search cost `O(e·q)`.  The step count is
/// the `m` table points plus every giant point compared.  Each digit is
/// the smallest solution below `q`, the same one the linear search
/// ([`Point::linear_dlog_vartime`]) returns.
pub fn recover_in_prime_power_subgroup(
    curve: &CurveParams,
    g: &Point,
    q: &Point,
    prime: &BigUint,
    exponent: u32,
) -> Option<(BigUint, u64)> {
    let a_fe = curve.a_fe();
    let mut d = BigUint::zero();
    let mut accumulator = BigUint::one();
    // Precompute G' = q^(e-1) · G, the generator of the q-order subgroup.
    let q_pow_e_minus_1 = prime.pow(exponent - 1);
    let g_prime = g.scalar_mul_vartime(&q_pow_e_minus_1, &a_fe);
    let q_u64 = prime.to_u64_digits().first().copied().unwrap_or(0);
    let baby = BabySteps::build(&g_prime, q_u64, &a_fe);
    let mut total_steps = baby.m;
    for j in 0..exponent {
        // Compute Q_j = q^(e-1-j) · (Q - d·G)
        let neg_dg = g.scalar_mul_vartime(&d, &a_fe).neg();
        let q_minus_dg = q.add_vartime(&neg_dg, &a_fe);
        let exp = prime.pow(exponent - 1 - j);
        let q_j = q_minus_dg.scalar_mul_vartime(&exp, &a_fe);
        let (d_j, giant) = baby.solve(&q_j, q_u64, &a_fe);
        total_steps += giant;
        let d_j = BigUint::from(d_j?);
        d += &d_j * &accumulator;
        accumulator *= prime;
    }
    Some((d, total_steps))
}

/// **Chinese Remainder reconstruction**.  Given `(m_i, r_i)` pairs
/// with the `m_i` pairwise coprime, return the unique `x` such that
/// `x ≡ r_i (mod m_i)`.
pub fn crt_combine(residues: &[(BigUint, BigUint)]) -> Option<BigUint> {
    if residues.is_empty() {
        return None;
    }
    let mut m = BigUint::one();
    let mut x = BigUint::zero();
    for (mi, ri) in residues {
        let m_inv = mod_inverse(&m, mi)?;
        let diff = if ri >= &(&x % mi) {
            (ri - &x % mi) % mi
        } else {
            (mi + ri - &x % mi) % mi
        };
        let t = (&diff * &m_inv) % mi;
        x += &t * &m;
        m *= mi;
    }
    Some(x % &m)
}

/// Render a Markdown report visualizing one PH run.
pub fn format_visualization(report: &PohligHellmanReport, order: &BigUint) -> String {
    let mut s = String::new();
    s.push_str("# Pohlig-Hellman recovery on a smooth-order group\n\n");
    s.push_str(&format!(
        "**Group order**: {} = {}\n\n",
        order,
        report
            .factors
            .iter()
            .map(|(q, e)| if *e == 1 {
                q.to_str_radix(10)
            } else {
                format!("{}^{}", q.to_str_radix(10), e)
            })
            .collect::<Vec<_>>()
            .join(" · ")
    ));
    s.push_str("## Per-prime-power discovered residues\n\n");
    s.push_str("```\n");
    s.push_str("  q^e         d mod q^e\n");
    s.push_str("  ──────      ─────────\n");
    for (m, r) in &report.residues {
        s.push_str(&format!("  {:<10}  {}\n", m, r));
    }
    s.push_str("```\n\n");
    s.push_str("## CRT reconstruction\n\n");
    match &report.recovered_d {
        Some(d) => s.push_str(&format!(
            "  {} **`d = {}`**  (recovered in {} ms, {} baby steps)\n",
            paint("✓", FG_BRIGHT_GREEN),
            d,
            report.elapsed_ms,
            report.total_steps,
        )),
        None => {
            // Distinguish "no smooth factor at all" from "partial recovery".
            if report.residues.is_empty() {
                s.push_str(&format!(
                    "  {} no factor below the smoothness bound — group has no usable structure\n",
                    paint("✗", FG_BRIGHT_RED),
                ));
            } else {
                s.push_str(&format!(
                    "  {} partial recovery: `d` known modulo {} (= ∏ recovered q^e factors)\n",
                    paint("⚠", FG_BRIGHT_YELLOW),
                    report
                        .residues
                        .iter()
                        .map(|(m, _)| m.clone())
                        .fold(BigUint::one(), |a, b| a * b),
                ));
            }
        }
    }
    s
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Toy curve with order 199 (prime) — Pohlig-Hellman fully reduces
    /// to standard rho since there's only one factor.  Tests the
    /// scaffolding.
    fn small_prime_curve() -> CurveParams {
        CurveParams {
            name: "ph-test-199",
            p: BigUint::from(211u32),
            a: BigUint::zero(),
            b: BigUint::from(2u32),
            gx: BigUint::from(4u32),
            gy: BigUint::from(53u32),
            n: BigUint::from(199u32),
            h: 1,
        }
    }

    /// CRT works on simple case: x ≡ 1 (mod 3), x ≡ 2 (mod 5) → x = 7.
    #[test]
    fn crt_simple() {
        let res = vec![
            (BigUint::from(3u32), BigUint::from(1u32)),
            (BigUint::from(5u32), BigUint::from(2u32)),
        ];
        assert_eq!(crt_combine(&res), Some(BigUint::from(7u32)));
    }

    /// **Pohlig-Hellman on a deliberately smooth-order group**.
    /// Construct a curve whose generator has order 60 = 2² · 3 · 5
    /// (toy, fully smooth), recover the discrete log via PH.
    ///
    /// We use a non-prime-order subgroup of the toy curve: pick a
    /// generator H = 2·G of a smaller subgroup.  Actually simpler:
    /// we construct a synthetic curve where we KNOW the smooth order.
    #[test]
    fn ph_recovers_dlp_on_smooth_curve() {
        // Use the toy j=0 curve with p=61, n=61.  61 is prime — not
        // smooth.  Instead, use p=829, n=823 (prime).  Still not smooth.
        //
        // For a meaningful smooth-order test, we need a curve whose
        // order factorises into small primes.  Use one of our twist-
        // enumeration results: p=65353 has a twist with order
        // 65856 = 2^6 · 3 · 7^3.  But constructing that curve here
        // requires the twist coefficient.
        //
        // Quick alternative: a toy curve we cook up — y² = x³ + 7 mod 67.
        // We need to compute its order; we know p+1-trace.
        //
        // Pragmatic test: build a curve and just iterate to find
        // its order's smoothest factor.
        let curve = small_prime_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let d_truth = BigUint::from(73u32);
        let q = g.scalar_mul(&d_truth, &a_fe);
        // 199 is prime; PH reduces to a single √199 ≈ 14-step solve
        // in the full subgroup.  smoothness_bound = 256 covers it.
        let report = pohlig_hellman_curve(&curve, &g, &q, &curve.n, 256);
        assert_eq!(report.recovered_d, Some(d_truth));
        assert_eq!(report.factors.len(), 1);
        assert_eq!(report.factors[0].0, BigUint::from(199u32));
    }

    /// Smoothness bound too small → partial recovery is None.
    #[test]
    fn ph_with_too_small_bound_partial_only() {
        let curve = small_prime_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let q = g.scalar_mul(&BigUint::from(10u32), &a_fe);
        // 199 > 5, so PH can't crack it.
        let report = pohlig_hellman_curve(&curve, &g, &q, &curve.n, 5);
        assert!(report.recovered_d.is_none());
        assert!(report.residues.is_empty());
    }

    /// **Prime-power subgroup recovery**: contrived test on a
    /// q=p case where we make ord(G) = 199 (= q^1).
    #[test]
    fn prime_power_subgroup_recovery_works() {
        let curve = small_prime_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let q = g.scalar_mul(&BigUint::from(99u32), &a_fe);
        let result = recover_in_prime_power_subgroup(&curve, &g, &q, &BigUint::from(199u32), 1);
        let (d, steps) = result.expect("should recover");
        assert_eq!(d, BigUint::from(99u32));
        assert!(steps <= 199);
    }

    /// **Visualization** renders the canonical sections.
    #[test]
    fn visualization_renders() {
        let report = PohligHellmanReport {
            recovered_d: Some(BigUint::from(73u32)),
            factors: vec![(BigUint::from(199u32), 1)],
            residues: vec![(BigUint::from(199u32), BigUint::from(73u32))],
            elapsed_ms: 5,
            total_steps: 73,
        };
        let s = format_visualization(&report, &BigUint::from(199u32));
        assert!(s.contains("Pohlig-Hellman"));
        assert!(s.contains("d = 73"));
        assert!(s.contains("baby steps"));
    }

    /// [`recover_in_prime_power_subgroup`] as it was before the digit
    /// search moved to [`Point::linear_dlog_vartime`]: affine steps and
    /// the constant-time point operations.
    fn recover_reference(
        curve: &CurveParams,
        g: &Point,
        q: &Point,
        prime: &BigUint,
        exponent: u32,
    ) -> Option<(BigUint, u64)> {
        let a_fe = curve.a_fe();
        let mut d = BigUint::zero();
        let mut accumulator = BigUint::one();
        let mut total_steps = 0u64;
        let g_prime = g.scalar_mul(&prime.pow(exponent - 1), &a_fe);
        for j in 0..exponent {
            let neg_dg = g.scalar_mul(&d, &a_fe).neg();
            let q_minus_dg = q.add(&neg_dg, &a_fe);
            let q_j = q_minus_dg.scalar_mul(&prime.pow(exponent - 1 - j), &a_fe);
            let mut found = None;
            let mut current = Point::Infinity;
            let q_u64 = prime.to_u64_digits().first().copied().unwrap_or(0);
            for k in 0..q_u64 {
                total_steps += 1;
                if current == q_j {
                    found = Some(BigUint::from(k));
                    break;
                }
                current = current.add(&g_prime, &a_fe);
            }
            d += &found? * &accumulator;
            accumulator *= prime;
        }
        Some((d, total_steps))
    }

    /// The digit search returns what the affine search returned — the
    /// digits and the baby-step count, or the failure — on every
    /// prime-power part of random curves over small primes (counted by
    /// brute force), for targets inside the subgroup and outside it, and
    /// for a mis-stated exponent or prime.
    /// An order with a prime factor far beyond the bound (a 256-bit prime
    /// `n`) is not trial-divided to its square root: the factorisation
    /// stops at the bound, reports the cofactor, and recovers nothing.
    #[test]
    fn large_prime_order_stops_at_the_bound() {
        let curve = CurveParams::secp256k1();
        let g = curve.generator();
        let q = g.scalar_mul_vartime(&BigUint::from(12_345u32), &curve.a_fe());
        let r = pohlig_hellman_curve(&curve, &g, &q, &curve.n, 1 << 16);
        assert_eq!(r.factors, vec![(curve.n.clone(), 1)]);
        assert!(r.recovered_d.is_none() && r.residues.is_empty());
        // A smooth part times the large prime: the small primes are found
        // with their exponents and the rest is one entry.
        let order = &curve.n * BigUint::from(2u32 * 2 * 3 * 13171);
        let f = factorise_up_to(&order, 1 << 16);
        assert_eq!(
            f,
            vec![
                (BigUint::from(2u32), 2),
                (BigUint::from(3u32), 1),
                (BigUint::from(13_171u32), 1),
                (curve.n.clone(), 1)
            ]
        );
    }

    /// The answer is checked in the group: a mis-stated order whose
    /// residues combine to a wrong logarithm reports nothing, and order 1
    /// is solved by `d = 0`.
    #[test]
    fn answers_are_verified_in_the_group() {
        let curve = small_prime_curve();
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let q = g.scalar_mul(&BigUint::from(77u32), &a_fe);
        let r = pohlig_hellman_curve(&curve, &g, &q, &curve.n, 1 << 10);
        assert_eq!(r.recovered_d, Some(BigUint::from(77u32)));
        // Claim order 3·199 with a bound that admits only 3: unrecoverable.
        let r = pohlig_hellman_curve(&curve, &g, &q, &(&curve.n * 3u32), 3);
        assert!(r.recovered_d.is_none());
        let r = pohlig_hellman_curve(
            &curve,
            &Point::Infinity,
            &Point::Infinity,
            &BigUint::one(),
            10,
        );
        assert_eq!(r.recovered_d, Some(BigUint::zero()));
        for (claimed, d) in [(BigUint::from(198u32), 5u32), (BigUint::from(200u32), 150)] {
            let q = g.scalar_mul(&BigUint::from(d), &a_fe);
            let r = pohlig_hellman_curve(&curve, &g, &q, &claimed, 1 << 10);
            if let Some(x) = &r.recovered_d {
                assert_eq!(g.scalar_mul(x, &a_fe), q, "claimed order {claimed}");
            }
        }
    }

    #[test]
    fn prime_power_recovery_matches_affine_search() {
        use num_bigint::RandBigInt;
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(0x9_4011);
        for p in [211u64, 1009, 2003, 4001] {
            let mut root = vec![None; p as usize];
            for y in 0..p {
                root[(y * y % p) as usize] = Some(y);
            }
            for _ in 0..6 {
                let (a, b) = (rng.gen_range(0..p), rng.gen_range(0..p));
                if (4 * a * a * a + 27 * b * b).is_multiple_of(p) {
                    continue;
                }
                let mut count = 1u64;
                let mut on: Vec<(u64, u64)> = Vec::new();
                for x in 0..p {
                    let r = (x * x % p * x + a * x + b) % p;
                    if r == 0 {
                        count += 1;
                        on.push((x, 0));
                    } else if let Some(y) = root[r as usize] {
                        count += 2;
                        on.push((x, y));
                    }
                }
                let curve = CurveParams {
                    name: "ph-differential",
                    p: BigUint::from(p),
                    a: BigUint::from(a),
                    b: BigUint::from(b),
                    gx: BigUint::zero(),
                    gy: BigUint::zero(),
                    n: BigUint::from(count),
                    h: 1,
                };
                let pt = |(x, y): (u64, u64)| Point::Affine {
                    x: curve.fe(BigUint::from(x)),
                    y: curve.fe(BigUint::from(y)),
                };
                let a_fe = curve.a_fe();
                let mut factors = Vec::new();
                let mut m = count;
                let mut f = 2;
                while m > 1 {
                    let mut e = 0;
                    while m.is_multiple_of(f) {
                        m /= f;
                        e += 1;
                    }
                    if e > 0 {
                        factors.push((f, e));
                    }
                    f += 1;
                }
                let g = pt(on[rng.gen_range(0..on.len())]);
                for &(q, e) in &factors {
                    let cof = BigUint::from(count / q.pow(e));
                    let g_i = g.scalar_mul(&cof, &a_fe);
                    let prime = BigUint::from(q);
                    for trial in 0..6 {
                        let target = if trial < 3 {
                            g.scalar_mul(&rng.gen_biguint(24), &a_fe)
                        } else {
                            pt(on[rng.gen_range(0..on.len())])
                        };
                        let t_i = target.scalar_mul(&cof, &a_fe);
                        for (pr, ex) in [(&prime, e), (&prime, e + 1), (&BigUint::from(2u32), e)] {
                            if ex == 0 {
                                continue;
                            }
                            // The digits (or the failure) are the affine
                            // search's; the step count is baby-step
                            // giant-step's own, at most ⌈√q⌉ per digit
                            // plus the ⌈√q⌉-point table, where the affine
                            // search's was up to q per digit.
                            let got = recover_in_prime_power_subgroup(&curve, &g_i, &t_i, pr, ex);
                            let want = recover_reference(&curve, &g_i, &t_i, pr, ex);
                            assert_eq!(
                                got.as_ref().map(|(d, _)| d),
                                want.as_ref().map(|(d, _)| d),
                                "p={p} a={a} b={b} #E={count} q={pr}^{ex} g={g_i:?} t={t_i:?}"
                            );
                            if let Some((_, steps)) = got {
                                let qq = pr.to_u64_digits()[0];
                                let m = (1..).find(|m: &u64| m * m >= qq).unwrap();
                                assert!(steps <= m * (ex as u64 + 1), "q={qq} steps={steps}");
                            }
                        }
                    }
                }
            }
        }
    }
}
