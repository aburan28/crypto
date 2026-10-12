//! Small curves on which the full-group solvers finish in milliseconds.
//!
//! The NIST curves are far too large for an unbounded rho, so the generic
//! machinery is exercised on these: two prime-order prime-field curves and
//! Koblitz curves `E_a / F_{2^m}` for small odd `m`, whose group orders come
//! from the Lucas recurrence on the trace of `E_a / F_2` and whose
//! prime-order subgroups are checked with Miller–Rabin at construction.

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::group::{BinaryGroup, PrimeGroup};
use super::koblitz_group_order;
use crate::binary_ecc::curve::{self as bcurve, BinaryCurve, BinaryPoint};
use crate::binary_ecc::f2m::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::curve_traits::arith::{factor, is_prime};
use crate::ecc::curve::CurveParams;

/// `y² = x³ + 3x + 6` over `F_10007`, generator `(0, 1973)` of prime order
/// `10039` (shared with `ecdlp_variants::demo_group_small`).
pub fn prime_small() -> PrimeGroup {
    PrimeGroup::new(CurveParams {
        name: "toy-p10039",
        p: BigUint::from(10_007u32),
        a: BigUint::from(3u32),
        b: BigUint::from(6u32),
        gx: BigUint::zero(),
        gy: BigUint::from(1973u32),
        n: BigUint::from(10_039u32),
        h: 1,
    })
}

/// `y² = x³ + 6x + 4` over `F_99013`, generator `(0, 2)` of prime order
/// `98893`.
pub fn prime_mid() -> PrimeGroup {
    PrimeGroup::new(CurveParams {
        name: "toy-p98893",
        p: BigUint::from(99_013u32),
        a: BigUint::from(6u32),
        b: BigUint::from(4u32),
        gx: BigUint::zero(),
        gy: BigUint::from(2u32),
        n: BigUint::from(98_893u32),
        h: 1,
    })
}

/// Irreducible trinomials / pentanomials for the toy extension degrees.
fn toy_irreducible(m: u32) -> Option<IrreduciblePoly> {
    let low: &[u32] = match m {
        17 => &[3, 0],
        19 => &[5, 2, 1, 0],
        23 => &[5, 0],
        29 => &[2, 0],
        31 => &[3, 0],
        37 => &[6, 4, 1, 0],
        41 => &[3, 0],
        47 => &[5, 0],
        _ => return None,
    };
    Some(IrreduciblePoly {
        degree: m,
        low_terms: low.to_vec(),
    })
}

/// The Koblitz curve `E_a : y² + xy = x³ + ax² + 1` over `F_{2^m}` for a
/// supported toy `m`, with a generator of the prime-order subgroup.
///
/// `Err` when `m` has no table entry or `#E_a / h` is not prime (the toy
/// table: `(17, 1)`, `(19, 0)`, `(19, 1)`, `(23, 0)`, `(23, 1)`, `(41, 0)`
/// have prime subgroup orders; `(17, 0)` has `#E/4 = 137 · 239` and is a
/// Pohlig–Hellman fixture via [`koblitz_full_group`]).
pub fn koblitz(m: u32, a: u8) -> Result<BinaryGroup, String> {
    let (curve, _n_full) = koblitz_curve(m, a, false)?;
    BinaryGroup::new(&format!("toy-k{m}a{a}"), curve)
}

/// The same curve with the generator's order set to the *full* group
/// order `#E_a(F_{2^m}) = h·n` and a generator of that order, so the
/// Pohlig–Hellman path has a composite order to work on.
pub fn koblitz_full_group(m: u32, a: u8) -> Result<BinaryGroup, String> {
    let (curve, _) = koblitz_curve(m, a, true)?;
    Ok(BinaryGroup::without_frobenius(
        &format!("toy-k{m}a{a}-full"),
        curve,
    ))
}

fn koblitz_curve(m: u32, a: u8, full_group: bool) -> Result<(BinaryCurve, BigUint), String> {
    let irr = toy_irreducible(m).ok_or_else(|| format!("no toy irreducible for m = {m}"))?;
    let n_full = koblitz_group_order(a, m);
    let h = BigUint::from(if a == 0 { 4u32 } else { 2u32 });
    let n = &n_full / &h;
    if !full_group && !is_prime(&n) {
        return Err(format!("#E_{a}(F_2^{m}) / {h} = {n} is not prime"));
    }
    let full_primes: Vec<BigUint> = if full_group {
        let f = factor(&n_full, 1 << 20);
        assert!(f.complete(), "toy order factors by trial division");
        f.primes.into_iter().map(|(q, _)| q).collect()
    } else {
        Vec::new()
    };
    let a_fe = if a == 0 {
        F2mElement::zero(m)
    } else {
        F2mElement::one(m)
    };
    let mut curve = BinaryCurve {
        m,
        irreducible: irr,
        a: a_fe,
        b: F2mElement::one(m),
        generator: BinaryPoint::Infinity,
        order: n.clone(),
        cofactor: h.clone(),
    };
    // Deterministic point search: x = 2, 3, 4, … lifted via the half-trace.
    let probe = BinaryGroup::without_frobenius("probe", curve.clone());
    let mut x_int = 2u64;
    loop {
        let x = F2mElement::from_biguint(&BigUint::from(x_int), m);
        x_int += 1;
        let Some(p) = probe.lift_x(&x) else {
            continue;
        };
        debug_assert!(curve.is_on_curve(&p));
        if full_group {
            // Need a point of order exactly h·n: the group is cyclic (one
            // F_2-rational 2-torsion point), so check [N/q]P ≠ O for each
            // prime q | N.
            let ok = full_primes
                .iter()
                .all(|q| bcurve::scalar_mul(&curve, &p, &(&n_full / q)) != BinaryPoint::Infinity);
            if ok {
                curve.generator = p;
                curve.order = n_full.clone();
                curve.cofactor = BigUint::one();
                return Ok((curve, n_full));
            }
        } else {
            let g = bcurve::scalar_mul(&curve, &p, &h);
            if g != BinaryPoint::Infinity {
                debug_assert_eq!(bcurve::scalar_mul(&curve, &g, &n), BinaryPoint::Infinity);
                curve.generator = g;
                return Ok((curve, n_full));
            }
        }
    }
}

/// Toy names the CLI and the dispatcher accept: `toy-p10039`,
/// `toy-p98893`, `toy-k<m>a<a>`, `toy-k<m>a<a>-full`.
pub fn toy_names() -> Vec<&'static str> {
    vec![
        "toy-p10039",
        "toy-p98893",
        "toy-k17a1",
        "toy-k19a0",
        "toy-k19a1",
        "toy-k23a0",
        "toy-k23a1",
        "toy-k41a0",
        "toy-k17a0-full",
        "toy-k19a0-full",
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecdlp_nist::group::EcdlpGroup;

    #[test]
    fn toy_koblitz_generators_have_the_declared_order() {
        for (m, a) in [(17u32, 1u8), (19, 0), (19, 1), (23, 0), (23, 1)] {
            let g = koblitz(m, a).unwrap_or_else(|e| panic!("toy-k{m}a{a}: {e}"));
            assert!(g.is_on_curve(&g.generator()));
            assert!(g.is_identity(&g.mul(&g.generator(), g.order())));
            assert!(g.frobenius.is_some(), "toy Koblitz carries Frobenius data");
        }
    }

    #[test]
    fn toy_full_group_generator_has_composite_order() {
        let g = koblitz_full_group(17, 0).expect("toy");
        let n = g.order().clone();
        assert_eq!(n, BigUint::from(130_972u32));
        assert!(g.is_identity(&g.mul(&g.generator(), &n)));
        assert!(!g.is_identity(&g.mul(&g.generator(), &(&n / 2u32))));
    }

    #[test]
    fn toy_table_rejects_composite_subgroup_orders() {
        assert!(koblitz(17, 0).is_err(), "#E_0/4 = 137·239 is composite");
    }
}
