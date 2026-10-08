//! **NIST P-224 hidden-weakness hypotheses: the checks this repository can run.**
//!
//! A ranked list of four hypotheses about how a 1999-era designer could
//! have weakened the NIST prime curves (a weak class recognisable from
//! `b`, the special prime, hidden endomorphism or isogeny structure, an
//! unpublished ECDLP algorithm) cites evidence from this repository.  Most
//! of that evidence was gathered on P-256.  This module computes the P-224
//! counterparts, and the cross-curve facts the ranking compares P-224 with:
//!
//! 1. the ANSI X9.62 / FIPS 186 seed-to-`c` derivation and the check
//!    `b²·c ≡ −27 (mod p)`, for all five NIST prime curves;
//! 2. the 2-adic valuation and the factorisation of `p − 1` and `p + 1`;
//! 3. the factorisation of `t² − 4p`, the conductor of `Z[π]` and the
//!    fundamental discriminant `D_K`;
//! 4. the factorisation of the quadratic twist's order;
//! 5. the embedding degree, exactly, when `n − 1` factors within budget.
//!
//! The residue test on `b` reuses [`b_seed_profile`]; [`null_coefficients`]
//! supplies its calibration panel.
//!
//! Factoring is trial division followed by Pollard–Brent under an
//! iteration budget, so a cofactor that does not split is reported as
//! unfactored, never as prime.  A prime factor below `3.3·10²⁴` is proved
//! (trial division, or Miller–Rabin on the first 13 prime bases, which is
//! deterministic there); above it, it is a Miller–Rabin probable prime on
//! twenty bases, and [`Certainty`] says which.
//!
//! Report: `research/p224_hypothesis_ranking_20261008/`.

use crate::cryptanalysis::b_seed_profile;
use crate::cryptanalysis::isogeny_walk::field::is_probable_prime;
use crate::ecc::curve::CurveParams;
use crate::hash::sha1::sha1;
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use serde::Serialize;
use std::collections::BTreeMap;

/// A NIST prime curve with the seed FIPS 186-4 Appendix D publishes for it.
pub struct SeededCurve {
    pub params: CurveParams,
    /// The 160-bit `SEED`, hex.
    pub seed_hex: &'static str,
    /// The `c` FIPS 186-4 prints next to the seed, hex.
    pub published_c_hex: &'static str,
}

/// P-192, P-224, P-256, P-384 and P-521 with their FIPS 186-4 seeds.
pub fn nist_seeded_curves() -> Vec<SeededCurve> {
    vec![
        SeededCurve {
            params: CurveParams::p192(),
            seed_hex: "3045ae6fc8422f64ed579528d38120eae12196d5",
            published_c_hex: "3099d2bbbfcb2538542dcd5fb078b6ef5f3d6fe2c745de65",
        },
        SeededCurve {
            params: CurveParams::p224(),
            seed_hex: "bd71344799d5c7fcdc45b59fa3b9ab8f6a948bc5",
            published_c_hex: "5b056c7e11dd68f40469ee7f3c7a7d74f7d121116506d031218291fb",
        },
        SeededCurve {
            params: CurveParams::p256(),
            seed_hex: "c49d360886e704936a6678e1139d26b7819f7e90",
            published_c_hex: "7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d",
        },
        SeededCurve {
            params: CurveParams::p384(),
            seed_hex: "a335926aa319a27a1d00896a6773a4827acdac73",
            published_c_hex: "79d1e655f868f02fff48dcdee14151ddb80643c1406d0ca10dfe6fc52009540a\
                              495e8042ea5f744f6e184667cc722483",
        },
        SeededCurve {
            params: CurveParams::p521(),
            seed_hex: "d09e8800291cb85396cc6717393284aaa0da64ba",
            published_c_hex: "b48bfa5f420a34949539d2bdfc264eeeeb077688e44fbf0ad8f6d0edb37bd6b5\
                              33281000518e19f1b9ffbe0fe9ed8a3c2200b8f875e523868c70c1e5bf55bad637",
        },
    ]
}

/// ANSI X9.62 A.3.3.1 (FIPS 186-4 D.6): the integer `c` a seed generates
/// for a prime field of `t = ⌈log₂ p⌉` bits.
///
/// `s = ⌊(t − 1)/160⌋`, `h = t − 160·s`.  `W₀` is the `h` rightmost bits
/// of `SHA-1(seed)` with its leftmost bit cleared; `Wᵢ = SHA-1(seed + i
/// mod 2^g)` for `i = 1 … s`, with `g` the seed's bit length; `c` is the
/// integer `W₀ ‖ W₁ ‖ … ‖ W_s`.
pub fn x962_c_from_seed(seed: &[u8], t: u64) -> BigUint {
    let s = (t - 1) / 160;
    let h = t - 160 * s;
    let w0 = BigUint::from_bytes_be(&sha1(seed));
    let mut c = w0 & ((BigUint::one() << (h - 1)) - 1u8);
    let g = seed.len() * 8;
    let z = BigUint::from_bytes_be(seed);
    let modulus = BigUint::one() << g;
    for i in 1..=s {
        let si = (&z + i) % &modulus;
        let mut bytes = si.to_bytes_be();
        let mut padded = vec![0u8; seed.len() - bytes.len()];
        padded.append(&mut bytes);
        c = (c << 160) | BigUint::from_bytes_be(&sha1(&padded));
    }
    c
}

/// The seed check for one curve.
#[derive(Clone, Debug, Serialize)]
pub struct SeedCheck {
    pub curve: String,
    pub seed_hex: String,
    pub c_hex: String,
    /// The derived `c` equals the one FIPS 186-4 prints.
    pub c_matches_published: bool,
    /// `b²·c ≡ −27 (mod p)`: `b` is a square root of `−27/c`.
    pub b2c_is_minus_27: bool,
}

/// Run [`x962_c_from_seed`] on a curve and test `b²·c ≡ −27 (mod p)`.
pub fn check_seed(curve: &SeededCurve) -> SeedCheck {
    let p = &curve.params.p;
    let seed = hex::decode(curve.seed_hex).expect("seed hex");
    let c = x962_c_from_seed(&seed, p.bits());
    let published = BigUint::parse_bytes(curve.published_c_hex.as_bytes(), 16).expect("c hex");
    let b2c = (&curve.params.b * &curve.params.b % p) * &c % p;
    SeedCheck {
        curve: curve.params.name.to_string(),
        seed_hex: curve.seed_hex.to_string(),
        c_hex: format!("{:x}", c),
        c_matches_published: c == published,
        b2c_is_minus_27: b2c == p - 27u8,
    }
}

// ── Integer factoring ───────────────────────────────────────────────────────

/// How a prime factor's primality is known.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Certainty {
    /// Proved: below the square of the trial bound with no smaller
    /// factor, or below `3.3·10²⁴`, where Miller–Rabin on the first 13
    /// prime bases is deterministic.
    Proved,
    /// Miller–Rabin on twenty bases only.
    Probable,
}

/// A factorisation found within a budget.
#[derive(Clone, Debug, Default)]
pub struct Factored {
    /// Prime factors with exponents.
    pub primes: BTreeMap<BigUint, (u32, Option<Certainty>)>,
    /// Composite cofactors that did not split within the budget.
    pub composite: Vec<BigUint>,
    /// Pollard–Brent iterations spent.
    pub rho_iterations: u64,
}

impl Factored {
    fn add(&mut self, q: BigUint, e: u32, how: Certainty) {
        let entry = self.primes.entry(q).or_insert((0, Some(how)));
        entry.0 += e;
    }

    /// Every factor is prime: nothing is left unfactored.
    pub fn complete(&self) -> bool {
        self.composite.is_empty()
    }

    /// The largest prime factor.
    pub fn largest_prime(&self) -> Option<&BigUint> {
        self.primes.keys().next_back()
    }

    /// The JSON form.
    pub fn report(&self) -> FactorReport {
        FactorReport {
            factors: self
                .primes
                .iter()
                .map(|(q, &(e, how))| PrimePower {
                    prime: q.to_string(),
                    exponent: e,
                    bits: bits_f64(q),
                    certainty: how.unwrap_or(Certainty::Probable),
                })
                .collect(),
            unfactored: self
                .composite
                .iter()
                .map(|c| Cofactor {
                    value: c.to_string(),
                    bits: bits_f64(c),
                })
                .collect(),
            complete: self.complete(),
            rho_iterations: self.rho_iterations,
        }
    }
}

/// One prime power in a [`FactorReport`].
#[derive(Clone, Debug, Serialize)]
pub struct PrimePower {
    pub prime: String,
    pub exponent: u32,
    /// `log₂` of the prime.
    pub bits: f64,
    pub certainty: Certainty,
}

/// A composite cofactor left unsplit.
#[derive(Clone, Debug, Serialize)]
pub struct Cofactor {
    pub value: String,
    /// `log₂` of the cofactor.
    pub bits: f64,
}

/// The serialisable form of [`Factored`].
#[derive(Clone, Debug, Serialize)]
pub struct FactorReport {
    pub factors: Vec<PrimePower>,
    pub unfactored: Vec<Cofactor>,
    pub complete: bool,
    pub rho_iterations: u64,
}

/// `log₂ x`, from the top 64 bits.
pub fn bits_f64(x: &BigUint) -> f64 {
    let b = x.bits();
    if b <= 64 {
        return (x.iter_u64_digits().next().unwrap_or(0) as f64).log2();
    }
    let top = (x >> (b - 64)).iter_u64_digits().next().unwrap_or(0) as f64;
    top.log2() + (b - 64) as f64
}

/// The bound below which Miller–Rabin on the first 13 prime bases is
/// deterministic (Sorenson–Webster 2015: `ψ₁₃ > 3.317·10²⁴`).
fn mr_deterministic_bound() -> BigUint {
    BigUint::parse_bytes(b"3317044064679887385961981", 10).unwrap()
}

/// Factor `n` by trial division below `trial_bound`, then Pollard–Brent
/// with at most `rho_budget` iterations in total.
pub fn factor(n: &BigUint, trial_bound: u32, rho_budget: u64) -> Factored {
    let mut out = Factored::default();
    if n.is_zero() {
        out.composite.push(BigUint::zero());
        return out;
    }
    let mut rest = n.clone();
    for q in b_seed_profile::primes_up_to(trial_bound as u64) {
        let qb = BigUint::from(q);
        if &qb * &qb > rest {
            break;
        }
        let mut e = 0u32;
        while (&rest % &qb).is_zero() {
            rest /= &qb;
            e += 1;
        }
        if e > 0 {
            out.add(qb, e, Certainty::Proved);
        }
    }
    let trial_sq = BigUint::from(trial_bound) * BigUint::from(trial_bound);
    let mr_bound = mr_deterministic_bound();
    let mut stack = vec![rest];
    while let Some(m) = stack.pop() {
        if m.is_one() {
            continue;
        }
        if m < trial_sq {
            out.add(m, 1, Certainty::Proved);
            continue;
        }
        if is_probable_prime(&m) {
            let how = if m < mr_bound {
                Certainty::Proved
            } else {
                Certainty::Probable
            };
            out.add(m, 1, how);
            continue;
        }
        let mut split = None;
        for c in 1..=8u64 {
            if out.rho_iterations >= rho_budget {
                break;
            }
            split = pollard_brent(&m, c, rho_budget, &mut out.rho_iterations);
            if split.is_some() {
                break;
            }
        }
        match split {
            Some(d) => {
                let other = &m / &d;
                stack.push(d);
                stack.push(other);
            }
            None => out.composite.push(m),
        }
    }
    out
}

/// Brent's variant of Pollard rho on `x ↦ x² + c mod n`, multiplying 128
/// differences per gcd.  Returns a proper factor, or `None` once `used`
/// reaches `budget` or the walk closes without one.
pub fn pollard_brent(n: &BigUint, c: u64, budget: u64, used: &mut u64) -> Option<BigUint> {
    if n.is_even() {
        return Some(BigUint::from(2u8));
    }
    let c = BigUint::from(c);
    let step = |x: &BigUint| (x * x + &c) % n;
    let one = BigUint::one();
    let block: u64 = 128;
    let mut y = BigUint::from(2u8);
    let mut x = y.clone();
    let mut ys = y.clone();
    let mut q = one.clone();
    let mut g = one.clone();
    let mut r: u64 = 1;
    while g.is_one() {
        if *used >= budget {
            return None;
        }
        x = y.clone();
        for _ in 0..r {
            y = step(&y);
        }
        *used += r;
        let mut k = 0;
        while k < r && g.is_one() {
            ys = y.clone();
            let m = block.min(r - k);
            for _ in 0..m {
                y = step(&y);
                let d = if x > y { &x - &y } else { &y - &x };
                q = q * d % n;
            }
            *used += m;
            g = q.gcd(n);
            k += m;
        }
        r *= 2;
    }
    if &g == n {
        // The block overshot: replay it one step at a time.
        g = one.clone();
        for _ in 0..=block {
            ys = step(&ys);
            let d = if x > ys { &x - &ys } else { &ys - &x };
            g = d.gcd(n);
            if !g.is_one() {
                break;
            }
        }
    }
    if g.is_one() || &g == n {
        None
    } else {
        Some(g)
    }
}

/// The 2-adic valuation of a nonzero `n`.
pub fn v2(n: &BigUint) -> u64 {
    n.trailing_zeros().unwrap_or(0)
}

/// The multiplicative order of `a` modulo `m`, given the complete
/// factorisation of `λ`, a multiple of that order.
pub fn multiplicative_order(a: &BigUint, m: &BigUint, lambda: &BigUint, f: &Factored) -> BigUint {
    assert!(f.complete(), "the order needs a complete factorisation");
    assert!(
        a.modpow(lambda, m).is_one(),
        "λ is not a multiple of the order"
    );
    let mut k = lambda.clone();
    for (q, &(e, _)) in &f.primes {
        for _ in 0..e {
            let t = &k / q;
            if a.modpow(&t, m).is_one() {
                k = t;
            } else {
                break;
            }
        }
    }
    k
}

// ── p ± 1 ───────────────────────────────────────────────────────────────────

/// The structure of `p − 1` or `p + 1` for one curve.
#[derive(Clone, Debug, Serialize)]
pub struct PrimeSide {
    pub curve: String,
    /// `"p-1"` or `"p+1"`.
    pub side: String,
    /// `v₂(p ∓ 1)`: the depth of a 2-power Kummer tower (`p − 1`) or of
    /// the split-torus Dickson tower (`p + 1`) over `F_p`.
    pub v2: u64,
    /// `log₂` of the largest prime factor found.
    pub largest_prime_bits: Option<f64>,
    pub factorization: FactorReport,
}

/// [`PrimeSide`] for both sides of one field prime.
pub fn prime_sides(name: &str, p: &BigUint, trial_bound: u32, rho_budget: u64) -> Vec<PrimeSide> {
    [("p-1", p - 1u8), ("p+1", p + 1u8)]
        .into_iter()
        .map(|(side, m)| {
            let f = factor(&m, trial_bound, rho_budget);
            PrimeSide {
                curve: name.to_string(),
                side: side.to_string(),
                v2: v2(&m),
                largest_prime_bits: f.largest_prime().map(bits_f64),
                factorization: f.report(),
            }
        })
        .collect()
}

// ── CM discriminant ─────────────────────────────────────────────────────────

/// `t² − 4p`, its conductor and fundamental discriminant.
#[derive(Clone, Debug, Serialize)]
pub struct CmDiscriminant {
    pub curve: String,
    /// The Frobenius trace `t = p + 1 − h·n`.
    pub trace: String,
    /// `Δ = t² − 4p`, negative.
    pub delta: String,
    pub delta_bits: f64,
    /// Factorisation of `|Δ|`.
    pub factorization: FactorReport,
    /// `f = [O_K : Z[π]]`, when `|Δ|` factored completely.
    pub conductor: Option<String>,
    /// `D_K = Δ / f²`, when `|Δ|` factored completely.
    pub fundamental: Option<String>,
    pub fundamental_bits: Option<f64>,
}

/// The Frobenius trace of a curve: `p + 1 − h·n`, as a signed integer.
pub fn trace(curve: &CurveParams) -> BigInt {
    let order = &curve.n * BigUint::from(curve.h);
    BigInt::from(&curve.p + 1u8) - BigInt::from(order)
}

/// Factor `t² − 4p` and split it into conductor² · fundamental discriminant.
pub fn cm_discriminant(curve: &CurveParams, trial_bound: u32, rho_budget: u64) -> CmDiscriminant {
    let t = trace(curve);
    let delta = &t * &t - BigInt::from(&curve.p * 4u8);
    let abs = delta.magnitude().clone();
    let f = factor(&abs, trial_bound, rho_budget);
    let (conductor, fundamental) = if f.complete() {
        // Δ ≡ t² (mod 4), so the odd part of the square is the conductor
        // when Δ ≡ 1 (mod 4); the curves here have odd t.
        let mut cond = BigUint::one();
        let mut squarefree = BigUint::one();
        for (q, &(e, _)) in &f.primes {
            cond *= q.pow(e / 2);
            if e % 2 == 1 {
                squarefree *= q;
            }
        }
        let mut dk = BigInt::from_biguint(Sign::Minus, squarefree);
        // −d with d squarefree is fundamental iff −d ≡ 1 (mod 4); otherwise
        // the field discriminant is −4d and 2 leaves the conductor.
        if dk.mod_floor(&BigInt::from(4)) != BigInt::one() {
            dk *= 4;
            cond /= 2u8;
        }
        (Some(cond), Some(dk))
    } else {
        (None, None)
    };
    CmDiscriminant {
        curve: curve.name.to_string(),
        trace: t.to_string(),
        delta: delta.to_string(),
        delta_bits: bits_f64(&abs),
        factorization: f.report(),
        conductor: conductor.as_ref().map(|c| c.to_string()),
        fundamental_bits: fundamental.as_ref().map(|d| bits_f64(d.magnitude())),
        fundamental: fundamental.map(|d| d.to_string()),
    }
}

/// Check a published factorisation: the factors multiply to `n` and each
/// is a (probable) prime.  Returns the certainty of the weakest factor,
/// or `None` if the product or a primality test fails.
pub fn verify_factorization(n: &BigUint, factors: &[&str]) -> Option<Certainty> {
    let mut product = BigUint::one();
    let mut weakest = Certainty::Proved;
    let bound = mr_deterministic_bound();
    for s in factors {
        let q = BigUint::parse_bytes(s.as_bytes(), 10)?;
        if !is_probable_prime(&q) {
            return None;
        }
        if q >= bound {
            weakest = Certainty::Probable;
        }
        product *= q;
    }
    (&product == n).then_some(weakest)
}

// ── Twist ───────────────────────────────────────────────────────────────────

/// The order of the quadratic twist and its factorisation.
#[derive(Clone, Debug, Serialize)]
pub struct TwistOrder {
    pub curve: String,
    /// `#E' = 2(p + 1) − #E = p + 1 + t`.
    pub order: String,
    pub factorization: FactorReport,
    /// `log₂` of the largest prime factor.
    pub largest_prime_bits: Option<f64>,
    /// `log₂ √(πℓ/4)`: expected rho cost in the largest prime-order
    /// subgroup `ℓ`, valid only when the factorisation is complete.
    pub rho_bits_largest: Option<f64>,
}

/// Factor the twist order `p + 1 + t`.
pub fn twist_order(curve: &CurveParams, trial_bound: u32, rho_budget: u64) -> TwistOrder {
    let order = &curve.p + 1u8 + &curve.p + 1u8 - &curve.n * BigUint::from(curve.h);
    let f = factor(&order, trial_bound, rho_budget);
    let largest = f.largest_prime().map(bits_f64);
    let rho = if f.complete() {
        largest.map(|b| 0.5 * (b + (std::f64::consts::PI / 4.0).log2()))
    } else {
        None
    };
    TwistOrder {
        curve: curve.name.to_string(),
        order: order.to_string(),
        factorization: f.report(),
        largest_prime_bits: largest,
        rho_bits_largest: rho,
    }
}

// ── Embedding degree ────────────────────────────────────────────────────────

/// The embedding degree `k`: the order of `p` in `(Z/n)^*`.
#[derive(Clone, Debug, Serialize)]
pub struct EmbeddingDegree {
    pub curve: String,
    /// Factorisation of `n − 1`.
    pub n_minus_1: FactorReport,
    /// `k`, when `n − 1` factored completely.
    pub k: Option<String>,
    /// `(n − 1)/k`, when `k` is known.
    pub index: Option<String>,
    pub k_bits: Option<f64>,
}

/// Compute the embedding degree exactly when `n − 1` factors in budget.
pub fn embedding_degree(curve: &CurveParams, trial_bound: u32, rho_budget: u64) -> EmbeddingDegree {
    let nm1 = &curve.n - 1u8;
    let f = factor(&nm1, trial_bound, rho_budget);
    let k = f
        .complete()
        .then(|| multiplicative_order(&(&curve.p % &curve.n), &curve.n, &nm1, &f));
    EmbeddingDegree {
        curve: curve.name.to_string(),
        n_minus_1: f.report(),
        index: k.as_ref().map(|k| (&nm1 / k).to_string()),
        k_bits: k.as_ref().map(bits_f64),
        k: k.map(|k| k.to_string()),
    }
}

// ── b-residue calibration panel ─────────────────────────────────────────────

/// `count` pseudo-random coefficients in `[0, p)`, for calibrating the
/// residue test: draw `i` is `SHA-1(label ‖ i ‖ 0) ‖ SHA-1(label ‖ i ‖ 1)`
/// read as a 320-bit integer and reduced mod `p`.
pub fn null_coefficients(label: &str, p: &BigUint, count: u32) -> Vec<BigUint> {
    (0..count)
        .map(|i| {
            let mut wide = Vec::with_capacity(40);
            for half in 0..2u8 {
                let mut msg = label.as_bytes().to_vec();
                msg.extend_from_slice(&i.to_be_bytes());
                msg.push(half);
                wide.extend_from_slice(&sha1(&msg));
            }
            BigUint::from_bytes_be(&wide) % p
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_nist_seed_generates_its_b() {
        for curve in nist_seeded_curves() {
            let check = check_seed(&curve);
            assert!(check.b2c_is_minus_27, "{}: b²c ≢ −27", check.curve);
            assert!(check.c_matches_published, "{}: c differs", check.curve);
        }
    }

    #[test]
    fn a_wrong_seed_does_not_generate_b() {
        let mut curve = nist_seeded_curves().remove(1);
        curve.seed_hex = "bd71344799d5c7fcdc45b59fa3b9ab8f6a948bc6";
        assert!(!check_seed(&curve).b2c_is_minus_27);
    }

    #[test]
    fn p224_p_minus_1_is_2_96_times_2_128_minus_1() {
        let p = CurveParams::p224().p;
        let m: BigUint = (BigUint::one() << 96u32) * ((BigUint::one() << 128u32) - 1u8);
        assert_eq!(&p - 1u8, m);
        let sides = prime_sides("P-224", &p, 1 << 20, 1 << 16);
        let minus = &sides[0];
        assert_eq!(minus.v2, 96);
        assert!(minus.factorization.complete);
        let primes: Vec<&str> = minus
            .factorization
            .factors
            .iter()
            .map(|f| f.prime.as_str())
            .collect();
        assert_eq!(
            primes,
            [
                "2",
                "3",
                "5",
                "17",
                "257",
                "641",
                "65537",
                "274177",
                "6700417",
                "67280421310721"
            ]
        );
        assert!(minus
            .factorization
            .factors
            .iter()
            .all(|f| f.certainty == Certainty::Proved));
    }

    #[test]
    fn two_adic_valuations_of_the_nist_primes() {
        let expect = [
            ("P-192", 1, 64),
            ("P-224", 96, 1),
            ("P-256", 1, 96),
            ("P-384", 1, 32),
            ("P-521", 1, 521),
        ];
        for (curve, (name, minus, plus)) in nist_seeded_curves().iter().zip(expect) {
            let p = &curve.params.p;
            assert_eq!(curve.params.name, name);
            assert_eq!(v2(&(p - 1u8)), minus, "{name} p−1");
            assert_eq!(v2(&(p + 1u8)), plus, "{name} p+1");
        }
    }

    #[test]
    fn p224_cm_discriminant_has_conductor_3() {
        let d = cm_discriminant(&CurveParams::p224(), 1 << 20, 1 << 20);
        assert_eq!(d.trace, "4733100108545601916421827343930821");
        assert!(d.factorization.complete);
        assert_eq!(d.conductor.as_deref(), Some("3"));
        // SafeCurves, "CM field discriminants".
        assert_eq!(
            d.fundamental.as_deref(),
            Some("-9493061114565352281698673660738078664961855212656825491744070162387")
        );
    }

    #[test]
    fn published_p256_discriminant_factorization_checks() {
        let c = CurveParams::p256();
        let t = trace(&c);
        let delta = (BigInt::from(&c.p * 4u8) - &t * &t).magnitude().clone();
        let factors = [
            "3",
            "5",
            "456597257999",
            "1428624589419343516204097",
            "46523541035814968339936406074986559003387",
        ];
        assert_eq!(
            verify_factorization(&delta, &factors),
            Some(Certainty::Probable)
        );
        assert_eq!(verify_factorization(&delta, &factors[1..]), None);
    }

    #[test]
    fn brent_splits_a_semiprime_and_order_is_exact() {
        let n = BigUint::from(1_000_003u64) * BigUint::from(998_244_353u64);
        let f = factor(&n, 1000, 1 << 20);
        assert!(f.complete());
        assert_eq!(f.primes.len(), 2);
        // ord of 2 mod 101 is 100.
        let m = BigUint::from(101u8);
        let lam = BigUint::from(100u8);
        let fl = factor(&lam, 100, 0);
        assert_eq!(
            multiplicative_order(&BigUint::from(2u8), &m, &lam, &fl),
            lam
        );
        assert_eq!(
            multiplicative_order(&BigUint::from(5u8), &m, &lam, &fl),
            BigUint::from(25u8)
        );
    }

    #[test]
    fn small_twist_factors_of_p224() {
        let t = twist_order(&CurveParams::p224(), 1 << 10, 0);
        let small: Vec<(&str, u32)> = t
            .factorization
            .factors
            .iter()
            .map(|f| (f.prime.as_str(), f.exponent))
            .collect();
        assert_eq!(small, [("3", 2), ("11", 1), ("47", 1)]);
        assert!(!t.factorization.complete);
        assert!(t.rho_bits_largest.is_none());
    }

    #[test]
    fn null_coefficients_are_distinct_and_reduced() {
        let p = CurveParams::p224().p;
        let draws = null_coefficients("test", &p, 16);
        assert_eq!(draws.len(), 16);
        assert!(draws.iter().all(|b| b < &p));
        let mut sorted = draws.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 16);
    }
}
