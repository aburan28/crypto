//! The number theory a trait record rests on: budgeted factoring that
//! keeps what it could not split, the split `t² − 4q = v²·d_K`, the
//! Kronecker symbol at a prime, multiplicative orders and Lucas
//! sequences.
//!
//! Every routine here is deterministic.  The factoring budget counts
//! iterations of Pollard's map, not seconds, so a rebuilt trait file is
//! byte-identical on any machine.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::Status;
use crate::cryptanalysis::isogeny_class_search::is_probable_prime;

/// Below this bound Miller–Rabin to the bases `2, 3, …, 41` is a proof
/// (ψ₁₃, Sorenson and Webster 2015); [`is_probable_prime`] runs every
/// one of those bases, so a prime it passes below the bound is proved.
pub const MR_PROOF_BOUND: u128 = 3_317_044_064_679_887_385_961_981;

/// Trial division runs to this bound before Pollard's rho starts.
const TRIAL_BOUND: u32 = 1 << 16;

/// Products of `|x − y|` accumulated between two gcds in Brent's rho.
const RHO_BATCH: u64 = 128;

/// Whether `n` passes Miller–Rabin to the bases `2, …, 53`.
pub fn is_prime(n: &BigUint) -> bool {
    is_probable_prime(&BigInt::from(n.clone()))
}

/// [`Status::Proved`] for a prime below [`MR_PROOF_BOUND`],
/// [`Status::Probable`] above it.
pub fn prime_status(p: &BigUint) -> Status {
    match p.to_u128() {
        Some(v) if v < MR_PROOF_BOUND => Status::Proved,
        _ => Status::Probable,
    }
}

/// The weaker of two statuses: a value derived from two others is only
/// as good as the worse of them.
pub fn weaker(a: Status, b: Status) -> Status {
    a.max(b)
}

/// A factorisation that keeps what the budget could not split.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Factored {
    /// Prime factors (proved or probable, see [`prime_status`]) with
    /// their exponents, ascending.
    pub primes: Vec<(BigUint, u32)>,
    /// Composites rho could not split, with their multiplicities: a
    /// composite found as a perfect power keeps its exponent, so
    /// `t² − 4q = 7·c²` still yields `d_K` with `c` unfactored.  Empty
    /// when the factorisation is complete.
    pub composites: Vec<(BigUint, u32)>,
}

impl Factored {
    pub fn complete(&self) -> bool {
        self.composites.is_empty()
    }

    /// The product of [`Self::composites`]: `1` when complete.
    pub fn residue(&self) -> BigUint {
        self.composites
            .iter()
            .fold(BigUint::one(), |acc, (c, e)| acc * c.pow(*e))
    }

    /// [`Status::Proved`] or [`Status::Probable`] for a complete
    /// factorisation, by the weakest prime in it; [`Status::Bounded`] for
    /// an incomplete one.
    pub fn status(&self) -> Status {
        if !self.complete() {
            return Status::Bounded;
        }
        self.prime_status()
    }

    /// The weakest primality status among the primes found.
    pub fn prime_status(&self) -> Status {
        self.primes
            .iter()
            .map(|(p, _)| prime_status(p))
            .fold(Status::Proved, weaker)
    }

    /// The largest prime factor, when the factorisation is complete.
    pub fn largest_prime(&self) -> Option<&BigUint> {
        if !self.complete() {
            return None;
        }
        self.primes.last().map(|(p, _)| p)
    }
}

fn small_primes() -> &'static [u32] {
    static PRIMES: OnceLock<Vec<u32>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let n = TRIAL_BOUND as usize;
        let mut sieve = vec![true; n + 1];
        sieve[0] = false;
        sieve[1] = false;
        let mut i = 2;
        while i * i <= n {
            if sieve[i] {
                let mut j = i * i;
                while j <= n {
                    sieve[j] = false;
                    j += i;
                }
            }
            i += 1;
        }
        (0..=n).filter(|&i| sieve[i]).map(|i| i as u32).collect()
    })
}

/// `n = root^k` for the smallest prime `k ≤ 7` that fits, if any.
fn perfect_power(n: &BigUint) -> Option<(BigUint, u32)> {
    for k in [2u32, 3, 5, 7] {
        let root = n.nth_root(k);
        if root > BigUint::one() && root.pow(k) == *n {
            return Some((root, k));
        }
    }
    None
}

/// Brent's variant of Pollard's rho on an odd composite `n` that is not
/// a perfect power, with `|x − y|` batched between gcds.  `budget` caps
/// the iterations of `x ↦ x² + c` per constant `c`; `None` means no
/// factor was found within it.
fn rho(n: &BigUint, budget: u64) -> Option<BigUint> {
    let f = |x: &BigUint, c: &BigUint| (x * x + c) % n;
    let diff = |a: &BigUint, b: &BigUint| if a > b { a - b } else { b - a };
    for c in 1u32..=4 {
        let c = BigUint::from(c);
        let mut y = BigUint::from(2u8);
        let mut x = y.clone();
        let mut ys = y.clone();
        let mut q = BigUint::one();
        let mut g = BigUint::one();
        let mut r: u64 = 1;
        let mut steps: u64 = 0;
        while g.is_one() {
            x = y.clone();
            for _ in 0..r {
                y = f(&y, &c);
            }
            let mut k = 0;
            while k < r && g.is_one() {
                ys = y.clone();
                let batch = RHO_BATCH.min(r - k);
                for _ in 0..batch {
                    y = f(&y, &c);
                    q = (q * diff(&x, &y)) % n;
                }
                g = q.gcd(n);
                k += batch;
            }
            steps += 2 * r;
            r *= 2;
            if g.is_one() && steps > budget {
                return None;
            }
        }
        if g == *n {
            // The batch overshot: replay it one step at a time.
            loop {
                ys = f(&ys, &c);
                g = diff(&x, &ys).gcd(n);
                if !g.is_one() {
                    break;
                }
            }
        }
        if g != *n {
            return Some(g);
        }
    }
    None
}

/// Factor `n ≥ 1`: trial division to 2¹⁶, then perfect powers and
/// Brent's rho, each composite getting `budget` iterations.  A composite
/// rho cannot split is kept in [`Factored::composites`], never reported
/// as a prime.
pub fn factor(n: &BigUint, budget: u64) -> Factored {
    assert!(!n.is_zero(), "factor(0)");
    let mut primes: BTreeMap<BigUint, u32> = BTreeMap::new();
    let mut rest = n.clone();
    for &p in small_primes() {
        let bp = BigUint::from(p);
        if &bp * &bp > rest {
            break;
        }
        let mut e = 0;
        while (&rest % &bp).is_zero() {
            rest /= &bp;
            e += 1;
        }
        if e > 0 {
            primes.insert(bp, e);
        }
    }
    let mut composites: BTreeMap<BigUint, u32> = BTreeMap::new();
    let mut stack = vec![(rest, 1u32)];
    while let Some((m, mult)) = stack.pop() {
        if m.is_one() {
            continue;
        }
        if is_prime(&m) {
            *primes.entry(m).or_default() += mult;
        } else if let Some((root, k)) = perfect_power(&m) {
            stack.push((root, mult * k));
        } else if let Some(d) = rho(&m, budget) {
            let other = &m / &d;
            stack.push((d, mult));
            stack.push((other, mult));
        } else {
            *composites.entry(m).or_default() += mult;
        }
    }
    Factored {
        primes: primes.into_iter().collect(),
        composites: composites.into_iter().collect(),
    }
}

/// `t² − 4q = v²·d_K` with `d_K` a fundamental discriminant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DiscSplit {
    /// The fundamental discriminant `d_K` of `Q(π)`; `None` when a
    /// composite of odd multiplicity stayed unfactored, whose square part
    /// is then unknown.
    pub cm_disc: Option<BigInt>,
    /// The conductor `v = [O_K : Z[π]]`, or, when `cm_disc` is `None`,
    /// the part of it the factorisation exposed: a lower bound.
    pub conductor: BigUint,
    /// Primes of the conductor (or of its known part).
    pub conductor_primes: Vec<(BigUint, u32)>,
    /// Composites of the conductor rho could not split, with exponents.
    pub conductor_composites: Vec<(BigUint, u32)>,
    /// Of `cm_disc` and `conductor` as values; [`Status::Bounded`] when
    /// `cm_disc` is `None`.
    pub status: Status,
}

/// Split a negative discriminant `delta ≡ 0, 1 (mod 4)` from the
/// factorisation of `|delta|`.  The 2-adic part is exact whatever stays
/// unfactored: trial division removed every factor of 2, and an odd
/// square is `1 (mod 8)`, so an unknown square part cannot change the
/// class of `d_K` modulo 4.
pub fn split_discriminant(delta: &BigInt, fac: &Factored) -> DiscSplit {
    assert_eq!(delta.sign(), Sign::Minus, "an imaginary quadratic order");
    let mut square_free = BigUint::one();
    let mut conductor_primes: Vec<(BigUint, u32)> = Vec::new();
    for (p, e) in &fac.primes {
        if e % 2 == 1 {
            square_free *= p;
        }
        if e / 2 > 0 {
            conductor_primes.push((p.clone(), e / 2));
        }
    }
    let mut odd_residue = BigUint::one();
    let mut conductor_composites: Vec<(BigUint, u32)> = Vec::new();
    for (c, e) in &fac.composites {
        if e % 2 == 1 {
            odd_residue *= c;
        }
        if e / 2 > 0 {
            conductor_composites.push((c.clone(), e / 2));
        }
    }
    // s = −(square-free part).  The odd residue's own square part is
    // 1 (mod 8), so s mod 4 is known exactly.
    let s_mod4 = {
        let known = (&square_free * &odd_residue) % 4u8;
        (4 - known.to_u32().expect("< 4")) % 4
    };
    // d_K = s when s ≡ 1 (mod 4); otherwise d_K = 4s and v loses a 2.
    if s_mod4 != 1 {
        let two = BigUint::from(2u8);
        let slot = conductor_primes
            .iter_mut()
            .find(|(p, _)| *p == two)
            .expect("t² − 4q ≡ 0 (mod 4) carries a square factor of 2");
        slot.1 -= 1;
        if slot.1 == 0 {
            conductor_primes.retain(|(p, _)| *p != two);
        }
        square_free *= 4u8;
    }
    let conductor = conductor_primes
        .iter()
        .chain(&conductor_composites)
        .fold(BigUint::one(), |acc, (p, e)| acc * p.pow(*e));
    let exact = odd_residue.is_one();
    DiscSplit {
        cm_disc: exact.then(|| -BigInt::from(square_free)),
        conductor,
        conductor_primes,
        conductor_composites,
        status: if exact {
            fac.prime_status()
        } else {
            Status::Bounded
        },
    }
}

/// The Kronecker symbol `(d / ℓ)` at a prime `ℓ`.
pub fn kronecker_prime(d: &BigInt, ell: u64) -> i32 {
    if ell == 2 {
        return match d.mod_floor(&BigInt::from(8u8)).to_u32().expect("< 8") {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        };
    }
    let a = d
        .mod_floor(&BigInt::from(ell))
        .to_u64()
        .expect("reduced mod ℓ");
    if a == 0 {
        return 0;
    }
    let e = BigUint::from(a).modpow(&BigUint::from((ell - 1) / 2), &BigUint::from(ell));
    if e.is_one() {
        1
    } else {
        -1
    }
}

/// What the Frobenius discriminant says at one small prime `ℓ`, exactly,
/// without factoring it: `ℓ`'s exponent in the conductor (the depth of
/// the `ℓ`-isogeny volcano) and how `ℓ` splits in `Q(π)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EllLocal {
    pub depth: u32,
    /// `(d_K / ℓ)`: `1` split, `−1` inert, `0` ramified.
    pub splitting: i32,
}

pub fn ell_local(delta: &BigInt, ell: u64) -> EllLocal {
    let bl = BigInt::from(ell);
    let mut e = 0u32;
    let mut u = delta.clone();
    while (&u % &bl).is_zero() {
        u /= &bl;
        e += 1;
    }
    if ell == 2 {
        // Δ = 2^e·u, u odd.  An even e with u ≡ 1 (mod 4) leaves d_K
        // odd; u ≡ 3 (mod 4) needs 4 | d_K; an odd e needs 8 | d_K.
        return if e % 2 == 0 && u.mod_floor(&BigInt::from(4u8)).is_one() {
            EllLocal {
                depth: e / 2,
                splitting: kronecker_prime(&u, 2),
            }
        } else if e % 2 == 0 {
            EllLocal {
                depth: e.saturating_sub(2) / 2,
                splitting: 0,
            }
        } else {
            EllLocal {
                depth: e.saturating_sub(3) / 2,
                splitting: 0,
            }
        };
    }
    let depth = e / 2;
    let rest = delta / bl.pow(2 * depth);
    EllLocal {
        depth,
        splitting: kronecker_prime(&rest, ell),
    }
}

/// The order of `a` modulo the prime `r`, from the factorisation of
/// `r − 1` (which must be complete).
pub fn multiplicative_order(a: &BigUint, r: &BigUint, r_minus_1: &Factored) -> BigUint {
    assert!(r_minus_1.complete());
    let a = a % r;
    let mut k = r - 1u8;
    for (p, e) in &r_minus_1.primes {
        for _ in 0..*e {
            let cand = &k / p;
            if a.modpow(&cand, r).is_one() {
                k = cand;
            } else {
                break;
            }
        }
    }
    k
}

/// The least `k ≤ limit` with `a^k ≡ 1 (mod r)`, if there is one.
pub fn order_at_most(a: &BigUint, r: &BigUint, limit: u64) -> Option<u64> {
    let a = a % r;
    let mut x = a.clone();
    for k in 1..=limit {
        if x.is_one() {
            return Some(k);
        }
        x = (x * &a) % r;
    }
    None
}

/// `V_m(s, q) = α^m + β^m` for the roots of `X² − sX + q`: the trace
/// over `GF(q^m)` of a curve of trace `s` over `GF(q)`.
pub fn lucas_v(s: &BigInt, q: &BigInt, m: u32) -> BigInt {
    let mut v0 = BigInt::from(2u8);
    let mut v1 = s.clone();
    if m == 0 {
        return v0;
    }
    for _ in 1..m {
        let v2 = s * &v1 - q * &v0;
        v0 = v1;
        v1 = v2;
    }
    v1
}

/// `log₂ |x|` as an `f64` for integers of any size.
pub fn log2_abs(x: &BigInt) -> f64 {
    let bits = x.bits();
    if bits <= 1000 {
        return x.abs().to_f64().expect("finite below 2^1000").log2();
    }
    let shift = bits - 64;
    (x.abs() >> shift).to_f64().expect("64 bits").log2() + shift as f64
}

/// Round to `digits` decimal places, so a written float is the same on
/// every platform that computes the same correctly rounded inputs.
pub fn round_to(x: f64, digits: i32) -> f64 {
    let s = 10f64.powi(digits);
    (x * s).round() / s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn big(s: &str) -> BigUint {
        s.parse().unwrap()
    }

    #[test]
    fn factor_is_complete_on_smooth_and_semiprime_inputs() {
        let n = big("1000000016000000063"); // 1000000007 · 1000000009
        let f = factor(&n, 1 << 20);
        assert!(f.complete());
        assert_eq!(
            f.primes,
            vec![(big("1000000007"), 1), (big("1000000009"), 1)]
        );
        assert_eq!(f.status(), Status::Proved);
        let f = factor(&BigUint::from(2u64.pow(10) * 3u64.pow(4) * 65537), 1);
        assert_eq!(
            f.primes,
            vec![
                (BigUint::from(2u8), 10),
                (BigUint::from(3u8), 4),
                (BigUint::from(65537u32), 1)
            ]
        );
    }

    #[test]
    fn factor_finds_large_squares_without_rho() {
        // ECC2K-130's conductor primes: |t² − 4q| = 7 · 263² · p², the
        // last factor far beyond rho's reach and found as a square root.
        let p = big("146505763881528721");
        let n = BigUint::from(7u8) * BigUint::from(263u32 * 263) * &p * &p;
        let f = factor(&n, 1);
        assert!(f.complete());
        assert_eq!(
            f.primes,
            vec![(BigUint::from(7u8), 1), (BigUint::from(263u32), 2), (p, 2)]
        );
    }

    #[test]
    fn an_exhausted_budget_leaves_a_residue_not_a_fake_prime() {
        let n = big("1000000016000000063");
        let f = factor(&n, 1);
        assert!(!f.complete());
        assert_eq!(f.composites, vec![(n.clone(), 1)]);
        assert_eq!(f.residue(), n);
        assert!(f.primes.is_empty());
        assert_eq!(f.status(), Status::Bounded);
    }

    #[test]
    fn an_unfactored_square_still_fixes_the_cm_discriminant() {
        // 7 · c² with c = 1000000007 · 1000000009 beyond a budget of 1.
        let c = big("1000000016000000063");
        let delta = -(BigInt::from(7u8) * BigInt::from(c.clone() * &c));
        let f = factor(delta.magnitude(), 1);
        assert_eq!(f.composites, vec![(c.clone(), 2)]);
        let s = split_discriminant(&delta, &f);
        assert_eq!(s.cm_disc, Some(BigInt::from(-7)));
        assert_eq!(s.conductor, c);
        assert_eq!(s.conductor_composites, vec![(c, 1)]);
        assert_eq!(s.status, Status::Proved);
        // An odd power leaves the square part unknown: d_K is not fixed.
        // (1000000009 · 1000000021 ≡ 1 (mod 4), so −7c is a discriminant.)
        let c1 = big("1000000030000000189");
        let delta = -(BigInt::from(7u8) * BigInt::from(c1));
        let s = split_discriminant(&delta, &factor(delta.magnitude(), 1));
        assert_eq!(s.cm_disc, None);
        assert_eq!(s.status, Status::Bounded);
    }

    #[test]
    fn discriminant_split_matches_hand_cases() {
        let case = |delta: i64| {
            let d = BigInt::from(delta);
            let f = factor(&d.magnitude().clone(), 1 << 16);
            let s = split_discriminant(&d, &f);
            (
                s.cm_disc.unwrap().to_i64().unwrap(),
                s.conductor.to_u64().unwrap(),
            )
        };
        assert_eq!(case(-7), (-7, 1));
        assert_eq!(case(-28), (-7, 2));
        assert_eq!(case(-4), (-4, 1));
        assert_eq!(case(-16), (-4, 2));
        assert_eq!(case(-8), (-8, 1));
        assert_eq!(case(-32), (-8, 2));
        assert_eq!(case(-12), (-3, 2));
        assert_eq!(case(-3 * 25 * 49), (-3, 35));
        assert_eq!(case(-20), (-20, 1));
        assert_eq!(case(-7 * 4 * 9), (-7, 6));
    }

    #[test]
    fn local_data_agrees_with_the_full_split() {
        for delta in (-4000i64..-2).filter(|d| d.rem_euclid(4) <= 1) {
            let d = BigInt::from(delta);
            let s = split_discriminant(&d, &factor(&d.magnitude().clone(), 1 << 16));
            let dk = s.cm_disc.unwrap();
            for ell in [2u64, 3, 5, 7, 11] {
                let loc = ell_local(&d, ell);
                let want = s
                    .conductor_primes
                    .iter()
                    .find(|(p, _)| *p == BigUint::from(ell))
                    .map_or(0, |(_, e)| *e);
                assert_eq!(loc.depth, want, "Δ = {delta}, ℓ = {ell}");
                assert_eq!(
                    loc.splitting,
                    kronecker_prime(&dk, ell),
                    "Δ = {delta}, ℓ = {ell}"
                );
            }
        }
    }

    #[test]
    fn multiplicative_order_and_lucas() {
        let r = BigUint::from(1_000_003u32);
        let f = factor(&BigUint::from(1_000_002u32), 1 << 16);
        let k = multiplicative_order(&BigUint::from(2u8), &r, &f);
        assert!(BigUint::from(2u8).modpow(&k, &r).is_one());
        assert_eq!(
            order_at_most(&BigUint::from(2u8), &r, 1_000_002).map(BigUint::from),
            Some(k)
        );
        // K_0 over GF(2): t₁ = −1; over GF(2^2), GF(2^3): −3, 5.
        let q = BigInt::from(2);
        assert_eq!(lucas_v(&BigInt::from(-1), &q, 2), BigInt::from(-3));
        assert_eq!(lucas_v(&BigInt::from(-1), &q, 3), BigInt::from(5));
    }
}
