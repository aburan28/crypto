//! Conductor-gap CM pairs: a crater curve and floor curves with the same
//! prime order and a certified prime conductor gap, with no isogeny computed.
//!
//! Fix `D = −19` (`h(D) = 1`, `j_D = −884736 = −96³`). For a prime `ℓ` and
//! an odd `t` with `q = (t² + 19ℓ²)/4` and `N = q + 1 − t` both prime,
//! every ordinary curve over `GF(q)` with trace `t` has `t² − 4q = −19ℓ²`,
//! so its endomorphism ring is an order of `Q(√−19)` whose conductor
//! divides `ℓ`: either `O_K` itself (the *crater*, the unique curve with
//! `j = j_D`, since the class number is one) or the order of conductor `ℓ`
//! (the *floor*, every other curve of order `N`). The two levels are joined
//! only by degree-`ℓ` isogenies.
//!
//! `N` is prime and exceeds the width of the Hasse interval, so one affine
//! point `P` with `[N]P = O` already proves `#E = N`; the constructor
//! returns that point as the generator. The same scalar multiplication is
//! the whole floor search: sample `(a, b)`, test `[N]P`, keep the curve or
//! its quadratic twist when the order is `N`, and reject `j = j_D`.
//!
//! Discriminants that cannot work, each by a congruence: `D = −7` (and any
//! `D ≡ 1 mod 8`) forces `t`, `f` even; `D = −11` forces `3 | qN`.
//! `research/conductor_gap_cm_pairs_20261010/` has the derivation and
//! Sage-verified pairs at `ℓ = 263, 8209, 131101`.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use crate::cryptanalysis::ic_boundary::{
    CountedGroup, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::residual_walk::is_prime_u64;

/// `|D|` for `D = −19`.
pub const ABS_D: u64 = 19;
/// `−j_D = 96³`, so `j_D ≡ q − 884736 (mod q)`.
pub const NEG_J_D: u64 = 884_736;

/// Which level of the volcano a curve sits on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GapLevel {
    /// `End(E) = O_K`, conductor 1: the one curve with `j = j_D`.
    Crater,
    /// `End(E)` of conductor `ℓ`: any other curve of order `N`.
    Floor,
}

impl std::fmt::Display for GapLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            GapLevel::Crater => "crater",
            GapLevel::Floor => "floor",
        })
    }
}

/// The field and group sizes a pair `(ℓ, t)` defines, once checked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GapParams {
    pub ell: u64,
    pub t: u64,
    pub q: u64,
    pub n: u64,
}

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}

fn submod(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        m - (b - a)
    }
}

fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}

/// Check `(ℓ, t)` and return `q` and `N`, or say which condition fails.
pub fn params(ell: u64, t: u64) -> Result<GapParams, String> {
    if ell < 3 || !is_prime_u64(ell) {
        return Err(format!("ell = {ell} must be an odd prime"));
    }
    if t.is_multiple_of(2) {
        return Err(format!("t = {t} must be odd and positive"));
    }
    let num = (t as u128) * (t as u128) + (ABS_D as u128) * (ell as u128) * (ell as u128);
    if !num.is_multiple_of(4) {
        return Err(format!("t^2 + 19 ell^2 = {num} is not divisible by 4"));
    }
    let q = num / 4;
    if q >= (1u128 << 62) {
        return Err(format!("q = {q} is not below 2^62"));
    }
    let q = q as u64;
    if !is_prime_u64(q) {
        return Err(format!("q = (t^2 + 19 ell^2)/4 = {q} is not prime"));
    }
    if t > q {
        return Err(format!("t = {t} is not below q + 1"));
    }
    let n = q + 1 - t;
    if !is_prime_u64(n) {
        return Err(format!("N = q + 1 - t = {n} is not prime"));
    }
    // One point of order N pins #E = N only when N exceeds the Hasse width.
    if (n as u128) * (n as u128) <= 16 * (q as u128) {
        return Err(format!(
            "N = {n} is not above 4 sqrt(q); the order certificate would not be unique"
        ));
    }
    Ok(GapParams { ell, t, q, n })
}

/// The crater curve: `j = j_D`, twisted so that its order is `N`, with a
/// point of order `N` as its generator.
pub fn crater(p: GapParams) -> Result<PrimeInstance, String> {
    let q = p.q;
    let j = submod(0, NEG_J_D % q, q);
    // y² = x³ + ax + b with a = 3j(1728 − j), b = 2j(1728 − j)².
    let k = submod(1728 % q, j, q);
    let a = mulmod(3, mulmod(j, k, q), q);
    let b = mulmod(2, mulmod(j, mulmod(k, k, q), q), q);
    let curve = PrimeCurve { p: q, a, b };
    let twist = quadratic_twist(&curve);
    if let Some(inst) = certify(curve, p, GapLevel::Crater) {
        return Ok(inst);
    }
    certify(twist, p, GapLevel::Crater).ok_or_else(|| {
        format!("neither j = -884736 nor its twist has order {} over GF({q}); parameters are inconsistent", p.n)
    })
}

/// A floor curve: a deterministic search from `seed` over `(a, b)`, each
/// candidate tested by one `[N]P`, the twist tried on a miss, `j = j_D`
/// rejected. About `19ℓ/4` samples per hit.
pub fn floor(p: GapParams, seed: u64) -> Result<PrimeInstance, String> {
    let q = p.q;
    let j_d = submod(0, NEG_J_D % q, q);
    let mut rng = StdRng::seed_from_u64(seed ^ (p.ell.rotate_left(32)) ^ p.t);
    // Expected 19ℓ/4 samples; allow a generous tail.
    let cap = 64 * ABS_D * p.ell + (1 << 16);
    for _ in 0..cap {
        let a = rng.gen_range(0..q);
        let b = rng.gen_range(1..q);
        let curve = PrimeCurve { p: q, a, b };
        let disc = addmod(
            mulmod(4, mulmod(mulmod(a, a, q), a, q), q),
            mulmod(27, mulmod(b, b, q), q),
            q,
        );
        if disc == 0 || j_invariant(&curve) == j_d {
            continue;
        }
        let twist = quadratic_twist(&curve);
        if let Some(inst) = certify(curve, p, GapLevel::Floor) {
            return Ok(inst);
        }
        if let Some(inst) = certify(twist, p, GapLevel::Floor) {
            return Ok(inst);
        }
    }
    Err(format!(
        "no floor curve found for ell = {}, t = {}, seed = {seed} within {cap} samples",
        p.ell, p.t
    ))
}

/// `j(E) = 1728 · 4a³ / (4a³ + 27b²)`.
pub fn j_invariant(curve: &PrimeCurve) -> u64 {
    let q = curve.p;
    let a3 = mulmod(4, mulmod(mulmod(curve.a, curve.a, q), curve.a, q), q);
    let den = addmod(a3, mulmod(27, mulmod(curve.b, curve.b, q), q), q);
    mulmod(mulmod(1728 % q, a3, q), powmod(den, q - 2, q), q)
}

/// The quadratic twist by the smallest non-residue `d`: `(a d², b d³)`.
fn quadratic_twist(curve: &PrimeCurve) -> PrimeCurve {
    let q = curve.p;
    let mut d = 2u64;
    while curve.legendre(d) != -1 {
        d += 1;
    }
    let d2 = mulmod(d, d, q);
    PrimeCurve {
        p: q,
        a: mulmod(curve.a, d2, q),
        b: mulmod(curve.b, mulmod(d2, d, q), q),
    }
}

/// The first affine point (smallest `x`, smaller `y`) with `[N]P = O`, if
/// the first point found has that order. Since `N` is prime and above the
/// Hasse width, one such point proves `#E = N`; one point without it proves
/// `#E ≠ N`, so a single test decides.
fn certify(curve: PrimeCurve, p: GapParams, level: GapLevel) -> Option<PrimeInstance> {
    let q = curve.p;
    let mut x = 0u64;
    let pt = loop {
        let rhs = curve.rhs(x);
        if rhs != 0 {
            if let Some(y) = curve.sqrt(rhs) {
                let y = y.min(q - y);
                break PrimePoint::affine(x, y);
            }
        }
        x += 1;
    };
    let mut scratch = GroupOps::default();
    if !curve.mul(&mut scratch, pt, p.n).infinity {
        return None;
    }
    debug_assert_eq!(
        j_invariant(&curve) == submod(0, NEG_J_D % q, q),
        level == GapLevel::Crater
    );
    let mut inst = PrimeInstance {
        name: String::new(),
        curve,
        group_order: p.n,
        r: p.n,
        cofactor: 1,
        generator: (pt.x, pt.y),
    };
    inst.name = inst.curve_id().slug;
    Some(inst)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ell_263_pair_matches_the_sage_verified_record() {
        // research/conductor_gap_cm_pairs_20261010/results/ell_263_q19.txt
        let p = params(263, 241).unwrap();
        assert_eq!((p.q, p.n), (343_073, 342_833));
        let c = crater(p).unwrap();
        assert_eq!(c.curve.p, 343_073);
        assert_eq!(j_invariant(&c.curve), 343_073 - 884_736 % 343_073);
        assert_eq!((c.group_order, c.r, c.cofactor), (p.n, p.n, 1));
        let f = floor(p, 7).unwrap();
        assert_ne!(j_invariant(&f.curve), j_invariant(&c.curve));
        assert_eq!(f.group_order, p.n);
        // Deterministic in the seed, different across seeds.
        let ab = |i: &PrimeInstance| (i.curve.a, i.curve.b, i.generator);
        assert_eq!(ab(&floor(p, 7).unwrap()), ab(&f));
        assert_ne!(ab(&floor(p, 8).unwrap()), ab(&f));
    }

    #[test]
    fn ell_8209_pair_builds() {
        let p = params(8209, 47).unwrap();
        assert_eq!((p.q, p.n), (320_092_037, 320_091_991));
        let c = crater(p).unwrap();
        assert_eq!(j_invariant(&c.curve), p.q - 884_736 % p.q);
    }

    #[test]
    fn bad_parameters_are_named() {
        assert!(params(263, 242).unwrap_err().contains("odd"));
        assert!(params(264, 241).unwrap_err().contains("prime"));
        // t = 3 at ell = 263: q = (9 + 19·263²)/4 = 328 583 = 7 · 46 941.
        assert!(params(263, 3).unwrap_err().contains("not prime"));
    }
}
