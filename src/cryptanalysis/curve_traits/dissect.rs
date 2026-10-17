//! The DiSSECT trait families (Sedlacek et al., "DiSSECT: Distinguisher
//! of Standard & Simulated Elliptic Curves via Traits") that follow from
//! the trace, the order and the endomorphism ring, under this module's
//! own definitions.  Parity with DiSSECT's code is not claimed: each
//! family below says what is computed, and `docs/curves/TRAITS.md` lists
//! the families left out and why.
//!
//! ## Frobenius mod `ℓ`
//!
//! Most families are properties of Frobenius acting on `E[ℓ] ≅ F_ℓ²`.
//! Its conjugacy class in `GL₂(F_ℓ)` is fixed by the characteristic
//! polynomial `X² − tX + q` unless that has a double root `λ`
//! (`ℓ | t² − 4q`); then it is the scalar `λ` when `ℓ | [End(E) : Z[π]]`
//! (the curve is above the floor of its `ℓ`-volcano) and a Jordan block
//! otherwise.  So with the volcano level in hand the class is exact, and
//! a representative matrix answers every question by enumeration:
//!
//! - `point_degree`: least `k` with a point of order `ℓ` over `GF(q^k)`;
//! - `full_degree`: least `k` with `E[ℓ] ⊆ E(GF(q^k))` (DiSSECT
//!   `torsion_extension`);
//! - `isogeny_degree`: least `k` with an `ℓ`-isogeny over `GF(q^k)`;
//!   `all_isogenies_degree`: least `k` with all `ℓ + 1` of them
//!   (`isogeny_extension`);
//! - `rational_isogenies`: `ℓ`-isogenies over `GF(q)`
//!   (`isogeny_neighbors`);
//! - `psi_degrees`: degrees of the irreducible factors of the division
//!   polynomial `ψ_ℓ` over `GF(q)`, as Frobenius orbits on
//!   `(E[ℓ] − O)/±1` (`division_polynomials`, `ℓ ≤ 5`).
//!
//! ## Factorisations
//!
//! `kn±1`, the traces `t_d` over `GF(q^d)`, and the twist orders over
//! `GF(q^d)` are factored by trial division to 2¹⁶, the committed factor
//! hints, 2¹⁰ rho iterations and a primality test of what is left (the
//! trace `t` itself gets the full budget).  Each is written
//! as `2^2·3·p250`: small primes in full, a large prime as `p<bits>`, an
//! unsplit composite as `c<bits>`.  A `c` means the factorisation is
//! partial.

use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

use super::arith::{self, lucas_v, Factored};
use super::{Factoring, Model, Status};

/// The primes `ℓ` the Frobenius-mod-`ℓ` families cover.
pub const TORSION_PRIMES: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19];

/// `ψ_ℓ` factor degrees are recorded for these.
pub const PSI_PRIMES: &[u64] = &[2, 3, 5];

/// Rho iterations for the families' factorisations, which are many and
/// mostly of numbers the size of `q`: enough to split off factors of
/// about 20 bits, no more.
const LIGHT_BUDGET: u64 = 1 << 10;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TorsionRecord {
    pub ell: u64,
    /// `scalar`, `jordan` (double root), `split`, `inert`.
    pub frobenius: Option<String>,
    pub point_degree: Option<u32>,
    pub full_degree: Option<u32>,
    pub isogeny_degree: Option<u32>,
    pub all_isogenies_degree: Option<u32>,
    pub rational_isogenies: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub psi_degrees: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DissectRecord {
    /// Of every exact value here: the order check's status (the trace
    /// everything rests on).  A torsion entry with no `frobenius` is
    /// unknown (a double root whose volcano level is unresolved).
    pub status: Status,
    pub torsion: Vec<TorsionRecord>,
    /// `|U_d(t, q)| = [Z[π] : Z[π^d]]`, so `v_d = v·|U_d|`: DiSSECT's
    /// conductor ratio for `d = 2, 3, 4`.
    pub conductor_ratio: BTreeMap<String, String>,
    /// `k·r ± 1` for `k = 1 … 8`, `r` the prime subgroup order.
    pub kn_factorization: BTreeMap<String, String>,
    /// `t_d = V_d(t, q)` for `d = 1, 2`.
    pub trace_factorization: BTreeMap<String, String>,
    /// `#E'(GF(q^d)) = q^d + 1 + t_d` for `d = 1, 2`.
    pub twist_order: BTreeMap<String, String>,
    /// `(r − 1) / ord_r(ℓ)` for `ℓ = 2 … 13`, when `r − 1` factors.
    pub small_prime_order_index: BTreeMap<String, String>,
    /// `r − 2^k` and `#E − 2^k` for the nearest power of two, as
    /// `2^k+<bits>` or `2^k-<bits>` (`2^k` when equal).
    pub pow_distance: BTreeMap<String, String>,
    /// `a³ / b² mod p` (prime fields, `b ≠ 0`): the X9.62 invariant.
    pub x962_invariant: Option<String>,
    /// For `w = 1, 2, 3` where at most [`HAMMING_MAX`] values are tried:
    /// how many `x` of Hamming weight `w` are x-coordinates of points.
    pub hamming_x: BTreeMap<String, u64>,
}

/// The most field elements `hamming_x` tries for one weight.
pub const HAMMING_MAX: u64 = 600;

// ── Frobenius mod ℓ ───────────────────────────────────────────────────

type Mat = [[u64; 2]; 2];

fn mat_mul(a: &Mat, b: &Mat, l: u64) -> Mat {
    let mut c = [[0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            c[i][j] = (a[i][0] * b[0][j] + a[i][1] * b[1][j]) % l;
        }
    }
    c
}

fn apply(m: &Mat, v: (u64, u64), l: u64) -> (u64, u64) {
    (
        (m[0][0] * v.0 + m[0][1] * v.1) % l,
        (m[1][0] * v.0 + m[1][1] * v.1) % l,
    )
}

fn is_scalar(m: &Mat) -> bool {
    m[0][1] == 0 && m[1][0] == 0 && m[0][0] == m[1][1]
}

fn vectors(l: u64) -> impl Iterator<Item = (u64, u64)> {
    (0..l)
        .flat_map(move |a| (0..l).map(move |b| (a, b)))
        .skip(1)
}

fn has_fixed(m: &Mat, l: u64) -> bool {
    vectors(l).any(|v| apply(m, v, l) == v)
}

/// `m·v` is a multiple of `v`.
fn eigen(m: &Mat, v: (u64, u64), l: u64) -> bool {
    let w = apply(m, v, l);
    (v.0 * w.1 % l + l - v.1 * w.0 % l).is_multiple_of(l)
}

fn has_eigenline(m: &Mat, l: u64) -> bool {
    vectors(l).any(|v| eigen(m, v, l))
}

fn eigenlines(m: &Mat, l: u64) -> u32 {
    let n = vectors(l).filter(|&v| eigen(m, v, l)).count() as u64;
    (n / (l - 1)) as u32
}

/// Orbit sizes of `m` on the non-zero vectors modulo `±1`, as
/// `1^2 2^5`.
fn psi_degrees(m: &Mat, l: u64) -> String {
    let canon = |v: (u64, u64)| {
        let n = ((l - v.0) % l, (l - v.1) % l);
        v.min(n)
    };
    let mut seen = std::collections::BTreeSet::new();
    let mut sizes: BTreeMap<u32, u32> = BTreeMap::new();
    for v in vectors(l) {
        let c = canon(v);
        if !seen.insert(c) {
            continue;
        }
        let mut size = 1;
        let mut w = canon(apply(m, c, l));
        while w != c {
            seen.insert(w);
            size += 1;
            w = canon(apply(m, w, l));
        }
        *sizes.entry(size).or_default() += 1;
    }
    sizes
        .iter()
        .map(|(d, n)| {
            if *n == 1 {
                d.to_string()
            } else {
                format!("{d}^{n}")
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn least_power(m: &Mat, l: u64, pred: impl Fn(&Mat) -> bool) -> u32 {
    let mut p = *m;
    // |GL₂(F_ℓ)| bounds every element's order.
    for k in 1..=(l * l * l * l) as u32 {
        if pred(&p) {
            return k;
        }
        p = mat_mul(&p, m, l);
    }
    unreachable!("every invertible matrix has finite order")
}

/// Frobenius mod `ℓ` for `ℓ` not the characteristic.  `scalar` decides
/// the double-root case: `Some(true)` above the floor, `Some(false)` on
/// it, `None` unknown.
pub fn torsion_record(
    ell: u64,
    trace: &BigInt,
    q: &BigUint,
    scalar: Option<bool>,
) -> TorsionRecord {
    let l = ell;
    let t = trace.mod_floor(&BigInt::from(l)).to_u64().expect("small");
    let qm = (q % l).to_u64().expect("small");
    let disc = (t * t + 4 * l * l - 4 * qm % l) % l;
    let double = if l == 2 { t == 0 } else { disc == 0 };
    let unknown = TorsionRecord {
        ell,
        frobenius: None,
        point_degree: None,
        full_degree: None,
        isogeny_degree: None,
        all_isogenies_degree: None,
        rational_isogenies: None,
        psi_degrees: None,
    };
    let (m, kind): (Mat, &str) = match (double, scalar) {
        (true, None) => return unknown,
        (true, Some(true)) => {
            let lambda = if l == 2 { 1 } else { t * l.div_ceil(2) % l };
            ([[lambda, 0], [0, lambda]], "scalar")
        }
        (true, Some(false)) => ([[0, (l - qm) % l], [1, t]], "jordan"),
        (false, _) => {
            // For ℓ = 2 a simple root means t odd, and X² + X + 1 is
            // irreducible over F₂.
            let kind = if l != 2 && arith::kronecker_prime(&BigInt::from(disc), l) == 1 {
                "split"
            } else {
                "inert"
            };
            ([[0, (l - qm) % l], [1, t]], kind)
        }
    };
    let ident = |p: &Mat| *p == [[1, 0], [0, 1]];
    TorsionRecord {
        ell,
        frobenius: Some(kind.into()),
        point_degree: Some(least_power(&m, l, |p| has_fixed(p, l))),
        full_degree: Some(least_power(&m, l, ident)),
        isogeny_degree: Some(least_power(&m, l, |p| has_eigenline(p, l))),
        all_isogenies_degree: Some(least_power(&m, l, is_scalar)),
        rational_isogenies: Some(eigenlines(&m, l)),
        psi_degrees: PSI_PRIMES.contains(&l).then(|| psi_degrees(&m, l)),
    }
}

// ── Factorisations ────────────────────────────────────────────────────

fn describe(n: &BigInt, fac: &Factored) -> String {
    if n.is_zero() {
        return "0".into();
    }
    let mut parts: Vec<String> = Vec::new();
    for (p, e) in &fac.primes {
        let base = if p.bits() <= 32 {
            p.to_string()
        } else {
            format!("p{}", p.bits())
        };
        parts.push(if *e == 1 { base } else { format!("{base}^{e}") });
    }
    for (c, e) in &fac.composites {
        let base = format!("c{}", c.bits());
        parts.push(if *e == 1 { base } else { format!("{base}^{e}") });
    }
    let body = if parts.is_empty() {
        "1".to_string()
    } else {
        parts.join("·")
    };
    if n.is_negative() {
        format!("-{body}")
    } else {
        body
    }
}

fn light(fx: &Factoring, n: &BigInt) -> String {
    if n.is_zero() {
        return "0".into();
    }
    describe(
        n,
        &arith::factor_hinted(n.magnitude(), LIGHT_BUDGET, &fx.hints),
    )
}

fn pow_distance(n: &BigUint) -> String {
    let k = n.bits() - 1;
    let lo = BigUint::one() << k;
    let hi = BigUint::one() << (k + 1);
    let (dl, dh) = (n - &lo, &hi - n);
    if dl.is_zero() {
        format!("2^{k}")
    } else if dl <= dh {
        format!("2^{k}+{}", dl.bits())
    } else {
        format!("2^{}-{}", k + 1, dh.bits())
    }
}

/// `U_d(t, q)`: `U_0 = 0`, `U_1 = 1`, `U_{d+1} = t·U_d − q·U_{d−1}`.
fn lucas_u(t: &BigInt, q: &BigInt, d: u32) -> BigInt {
    let (mut u0, mut u1) = (BigInt::zero(), BigInt::one());
    for _ in 1..d {
        let u2 = t * &u1 - q * &u0;
        u0 = u1;
        u1 = u2;
    }
    if d == 0 {
        u0
    } else {
        u1
    }
}

fn hamming_prime(p: &BigUint, a: &BigUint, b: &BigUint, w: u32) -> Option<u64> {
    let bits = p.bits();
    if binomial(bits, w)? > HAMMING_MAX {
        return None;
    }
    let half = (p - 1u8) >> 1;
    let mut count = 0;
    for_weight(bits, w, &mut |x: BigUint| {
        if &x >= p {
            return;
        }
        let f = (&x * &x % p * &x + a * &x + b) % p;
        if f.is_zero() || f.modpow(&half, p).is_one() {
            count += 1;
        }
    });
    Some(count)
}

fn hamming_binary(
    n: u32,
    irr: &crate::binary_ecc::IrreduciblePoly,
    a: &crate::binary_ecc::F2mElement,
    b: &crate::binary_ecc::F2mElement,
    w: u32,
) -> Option<u64> {
    use crate::binary_ecc::F2mElement;
    if binomial(u64::from(n), w)? > HAMMING_MAX {
        return None;
    }
    let mut count = 0;
    for_weight(u64::from(n), w, &mut |x: BigUint| {
        let x = F2mElement::from_biguint(&x, n);
        let inv = x.flt_inverse(irr).expect("weight ≥ 1");
        let s = x.add(a).add(&b.mul(&inv.square(irr), irr));
        if s.trace(irr) == 0 {
            count += 1;
        }
    });
    Some(count)
}

fn binomial(n: u64, k: u32) -> Option<u64> {
    let mut c: u64 = 1;
    for i in 0..u64::from(k) {
        c = c.checked_mul(n - i)? / (i + 1);
    }
    Some(c)
}

fn for_weight(bits: u64, w: u32, f: &mut dyn FnMut(BigUint)) {
    fn rec(start: u64, bits: u64, left: u32, acc: &BigUint, f: &mut dyn FnMut(BigUint)) {
        if left == 0 {
            f(acc.clone());
            return;
        }
        for i in start..bits {
            let mut next = acc.clone();
            next.set_bit(i, true);
            rec(i + 1, bits, left - 1, &next, f);
        }
    }
    rec(0, bits, w, &BigUint::zero(), f);
}

/// Every family for one curve.  `scalar(ℓ)` says whether Frobenius is a
/// scalar on `E[ℓ]` (see the module docs); `r` is the prime subgroup
/// order when known.
pub fn compute(
    fx: &Factoring,
    model: &Model,
    characteristic: &BigUint,
    trace: &BigInt,
    order: &BigUint,
    r: Option<&BigUint>,
    status: Status,
    scalar: &dyn Fn(u64) -> Option<bool>,
) -> DissectRecord {
    let q = model.q();
    let qi = BigInt::from(q.clone());
    let torsion = TORSION_PRIMES
        .iter()
        .filter(|&&l| *characteristic != BigUint::from(l))
        .map(|&l| torsion_record(l, trace, &q, scalar(l)))
        .collect();
    let conductor_ratio = (2..=4)
        .map(|d| (d.to_string(), lucas_u(trace, &qi, d).abs().to_string()))
        .collect();
    let mut kn_factorization = BTreeMap::new();
    let mut small_prime_order_index = BTreeMap::new();
    let mut pow = BTreeMap::new();
    pow.insert("order".to_string(), pow_distance(order));
    if let Some(r) = r {
        let ri = BigInt::from(r.clone());
        for k in 1..=8u32 {
            for (sign, delta) in [("+", 1i32), ("-", -1)] {
                let v = &ri * k + delta;
                kn_factorization.insert(format!("{k}n{sign}1"), light(fx, &v));
            }
        }
        pow.insert("subgroup".to_string(), pow_distance(r));
        if r > &BigUint::from(13u8) {
            let rm1 = fx.factor(&(r - 1u8));
            if rm1.complete() {
                for ell in [2u64, 3, 5, 7, 11, 13] {
                    let ord = arith::multiplicative_order(&BigUint::from(ell), r, &rm1);
                    small_prime_order_index.insert(ell.to_string(), ((r - 1u8) / ord).to_string());
                }
            }
        }
    }
    let mut trace_factorization = BTreeMap::new();
    let mut twist_order = BTreeMap::new();
    for d in 1..=2u32 {
        let td = lucas_v(trace, &qi, d);
        // t itself is no larger than 2√q: give it the full budget.
        let fac = if d == 1 && !td.is_zero() {
            describe(&td, &fx.factor(td.magnitude()))
        } else {
            light(fx, &td)
        };
        trace_factorization.insert(d.to_string(), fac);
        let twist = qi.pow(d) + 1 + &td;
        twist_order.insert(d.to_string(), light(fx, &twist));
    }
    let mut hamming_x = BTreeMap::new();
    let mut x962_invariant = None;
    for w in 1..=3u32 {
        let count = match model {
            Model::Prime { p, a, b } => hamming_prime(p, a, b, w),
            Model::Binary { n, irr, a, b, .. } => hamming_binary(*n, irr, a, b, w),
            Model::Other { .. } | Model::Extension { .. } => None,
        };
        if let Some(c) = count {
            hamming_x.insert(w.to_string(), c);
        }
    }
    if let Model::Prime { p, a, b } = model {
        if !b.is_zero() {
            let b2 = b * b % p;
            let inv = b2.modpow(&(p - 2u8), p);
            x962_invariant = Some((a * a % p * a % p * inv % p).to_string());
        }
    }
    DissectRecord {
        status,
        torsion,
        conductor_ratio,
        kn_factorization,
        trace_factorization,
        twist_order,
        small_prime_order_index,
        pow_distance: pow,
        x962_invariant,
        hamming_x,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frobenius_mod_ell_families_match_known_cases() {
        // t ≡ 2, q ≡ 1 mod 3 with Frobenius scalar: E[3] rational.
        let r = torsion_record(3, &BigInt::from(2), &BigUint::from(7u8), Some(true));
        assert_eq!(r.frobenius.as_deref(), Some("scalar"));
        assert_eq!(r.full_degree, Some(1));
        assert_eq!(r.rational_isogenies, Some(4));
        assert_eq!(r.psi_degrees.as_deref(), Some("1^4"));
        // The same class on the floor: one isogeny, full torsion only at k = 3.
        let j = torsion_record(3, &BigInt::from(2), &BigUint::from(7u8), Some(false));
        assert_eq!(j.rational_isogenies, Some(1));
        assert_eq!(j.point_degree, Some(1));
        assert_eq!(j.full_degree, Some(3));
        // Unknown level: nothing claimed.
        let u = torsion_record(3, &BigInt::from(2), &BigUint::from(7u8), None);
        assert_eq!(u.frobenius, None);
        // The Koblitz Frobenius over GF(2): X² + X + 2, irreducible mod 3.
        let k = torsion_record(3, &BigInt::from(-1), &BigUint::from(2u8), None);
        assert_eq!(k.frobenius.as_deref(), Some("inert"));
        assert_eq!(k.rational_isogenies, Some(0));
    }

    #[test]
    fn lucas_u_gives_the_index_of_z_pi_d() {
        // t = 1, q = 2: U_2 = 1, U_3 = −1, U_4 = −3.
        let (t, q) = (BigInt::one(), BigInt::from(2));
        assert_eq!(lucas_u(&t, &q, 2), BigInt::one());
        assert_eq!(lucas_u(&t, &q, 3), BigInt::from(-1));
        assert_eq!(lucas_u(&t, &q, 4), BigInt::from(-3));
        // v_d² d_K = t_d² − 4q^d with v_d = v·|U_d|: check d = 2, 3.
        for d in 2..=3 {
            let td = lucas_v(&t, &q, d);
            let lhs = &td * &td - BigInt::from(4) * q.pow(d);
            let u = lucas_u(&t, &q, d);
            assert_eq!(lhs, &u * &u * BigInt::from(-7));
        }
    }

    #[test]
    fn partial_factorisations_are_written_compactly() {
        let fx = Factoring::new(0, Vec::new());
        assert_eq!(light(&fx, &BigInt::from(360)), "2^3·3^2·5");
        let p = BigInt::from(1_000_000_007u64) * BigInt::from(1_000_000_009u64);
        assert_eq!(light(&fx, &(BigInt::from(12) * &p)), "2^2·3·c60");
        assert_eq!(pow_distance(&BigUint::from(1024u32)), "2^10");
        assert_eq!(pow_distance(&BigUint::from(1023u32)), "2^10-1");
    }
}
