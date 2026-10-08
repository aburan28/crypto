//! Checking the group order every other trait is derived from.
//!
//! The registry records `#E(GF(q))`; the trace, the Frobenius
//! discriminant and everything after them follow from it.  Two checks
//! certify it:
//!
//! - **Exhaustive count** for `q ≤ 2^COUNT_MAX_BITS`: every `x` is
//!   visited.  Proved.
//! - **Descent count** for a binary curve defined over a subfield
//!   `GF(2^k)`, `k ≤ subfield::COUNT_MAX_K`: count the descended model
//!   over `GF(2^k)`, then `t_n = V_{n/k}(t_k, 2^k)` (the Frobenius over
//!   `GF(2^n)` is the `n/k`-th power of the one over `GF(2^k)`).
//!   Proved.  This is how every Koblitz order is certified.
//! - **Generator certificate** from a recorded representation: `G` on
//!   the curve, `G ≠ O`, `[r]G = O` with `r` prime, so `r | #E`.  When
//!   `r > 4√q` the Hasse interval holds one multiple of `r`, which is
//!   then `#E`; the certificate is as strong as `r`'s primality.  When
//!   `r ≤ 4√q` it only shows the order is consistent.
//!
//! A count that disagrees, or a valid generator whose order is not `r`,
//! is an error in the registry, not a trait: the build stops.

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

use num_bigint::BigInt;

use super::arith::{is_prime, lucas_v, prime_status};
use super::{subfield, Model, Representation, Status};
use crate::binary_ecc::curve::{scalar_mul, BinaryCurve, BinaryPoint};
use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::semaev_decomp::Gf2;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;

/// Fields up to this many bits are counted point by point.
pub const COUNT_MAX_BITS: u64 = 22;

/// `#E` for `y² + xy = x³ + ax² + b` over `GF(2^n)`, `n ≤ 63`.
pub fn count_binary(gf: &Gf2, a: u64, b: u64) -> u64 {
    let n = gf.n;
    // Tr is F₂-linear: bit i of the mask is Tr(z^i).
    let mut tmask = 0u64;
    for i in 0..n {
        let mut x = 1u64 << i;
        let mut acc = 0u64;
        for _ in 0..n {
            acc ^= x;
            x = gf.sqr(x);
        }
        debug_assert!(acc <= 1);
        tmask |= acc << i;
    }
    let xs: Vec<u64> = (1..1u64 << n).collect();
    let mut inv = xs.clone();
    gf.batch_inv(&mut inv, &mut Vec::new());
    // O and (0, √b), then two points for each x ≠ 0 with
    // Tr(x + a + b/x²) = 0.
    let mut points = 2u64;
    for (x, xi) in xs.iter().zip(&inv) {
        let w = x ^ a ^ gf.mul(b, gf.sqr(*xi));
        if (w & tmask).count_ones() % 2 == 0 {
            points += 2;
        }
    }
    points
}

/// `#E` for `y² = x³ + ax + b` over `GF(p)`, `p < 2^32`.
pub fn count_prime(p: u64, a: u64, b: u64) -> u64 {
    let mut square = vec![false; p as usize];
    for y in 1..p {
        square[(y * y % p) as usize] = true;
    }
    let mut points = 1u64;
    for x in 0..p {
        let f = ((x * x % p) * x % p + a * x % p + b) % p;
        points += match f {
            0 => 1,
            _ if square[f as usize] => 2,
            _ => 0,
        };
    }
    points
}

/// The order check's outcome: a status and how it was reached.
pub type Outcome = (Status, &'static str);

/// Check `order` against the model.  `Err` when the model contradicts it.
pub fn check_order(
    model: &Model,
    order: &BigUint,
    reps: &[Representation],
) -> Result<Outcome, String> {
    let q = model.q();
    if model.size_bits() <= COUNT_MAX_BITS {
        let counted = match model {
            Model::Binary { irr, a, b, .. } => {
                let gf = Gf2::new(irr);
                count_binary(&gf, gf.from_element(a), gf.from_element(b))
            }
            Model::Prime { p, a, b } => count_prime(
                p.to_u64().expect("small"),
                a.to_u64().expect("reduced"),
                b.to_u64().expect("reduced"),
            ),
        };
        if BigUint::from(counted) != *order {
            return Err(format!(
                "exhaustive count gives #E = {counted}, the registry records {order}"
            ));
        }
        return Ok((Status::Proved, "exhaustive point count"));
    }
    if let Model::Binary { n, irr, a, b, .. } = model {
        let (_, k, tr_a) = subfield::degrees(*n, irr, a, b);
        if k < *n && k <= subfield::COUNT_MAX_K {
            let t_k = subfield::descended_trace_by_count(*n, k, irr, b, tr_a);
            let t_n = lucas_v(&t_k, &(BigInt::one() << k), n / k);
            let t_recorded = BigInt::from(q.clone()) + 1u8 - BigInt::from(order.clone());
            if t_n != t_recorded {
                return Err(format!(
                    "the model descends to GF(2^{k}) with t_k = {t_k}, so t = {t_n}; the registry records {t_recorded}"
                ));
            }
            return Ok((
                Status::Proved,
                "descent count: #E over GF(2^k), lifted by the Lucas sequence",
            ));
        }
    }
    let mut best: Option<Outcome> = None;
    for rep in reps {
        let Some(outcome) = generator_certificate(model, order, &q, rep)? else {
            continue;
        };
        if best.is_none_or(|b| outcome.0 < b.0) {
            best = Some(outcome);
        }
    }
    Ok(best.unwrap_or((
        Status::Unknown,
        "recorded order; no usable generator and the field is too large to count",
    )))
}

/// `None` when the representation's generator is not on this model (a
/// representation of another model cannot certify this one).
fn generator_certificate(
    model: &Model,
    order: &BigUint,
    q: &BigUint,
    rep: &Representation,
) -> Result<Option<Outcome>, String> {
    let r = &rep.subgroup_order;
    if &rep.cofactor * r != *order {
        return Err(format!(
            "cofactor {} · subgroup order {r} ≠ recorded order {order}",
            rep.cofactor
        ));
    }
    let (gx, gy) = &rep.generator;
    let killed = match model {
        Model::Binary { n, irr, a, b, .. } => {
            let g = BinaryPoint::Affine {
                x: F2mElement::from_biguint(gx, *n),
                y: F2mElement::from_biguint(gy, *n),
            };
            let curve = BinaryCurve {
                m: *n,
                irreducible: irr.clone(),
                a: a.clone(),
                b: b.clone(),
                generator: g.clone(),
                order: r.clone(),
                cofactor: rep.cofactor.clone(),
            };
            if gx.bits() > u64::from(*n) || gy.bits() > u64::from(*n) || !curve.is_on_curve(&g) {
                return Ok(None);
            }
            scalar_mul(&curve, &g, r) == BinaryPoint::Infinity
        }
        Model::Prime { p, a, b } => {
            if gx >= p || gy >= p {
                return Ok(None);
            }
            let lhs = gy * gy % p;
            let rhs = (gx * gx % p * gx + a * gx + b) % p;
            if lhs != rhs {
                return Ok(None);
            }
            let g = Point::Affine {
                x: FieldElement::new(gx.clone(), p.clone()),
                y: FieldElement::new(gy.clone(), p.clone()),
            };
            g.scalar_mul_vartime(r, &FieldElement::new(a.clone(), p.clone())) == Point::Infinity
        }
    };
    if !killed {
        return Err(format!("[r]G ≠ O for the recorded generator and r = {r}"));
    }
    if r.is_zero() || r.is_one() || !is_prime(r) {
        return Err(format!("recorded subgroup order {r} is not prime"));
    }
    // One multiple of r in [q + 1 − 2√q, q + 1 + 2√q] iff r > 4√q.
    if r * r > BigUint::from(16u8) * q {
        Ok(Some((
            prime_status(r),
            "generator certificate: r prime, [r]G = O, r > 4√q",
        )))
    } else {
        Ok(Some((
            Status::Bounded,
            "generator consistent: [r]G = O, but r ≤ 4√q leaves other multiples in the Hasse interval",
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::IrreduciblePoly;

    #[test]
    fn counts_match_known_orders() {
        // K_0 over GF(2^7), modulus z^7 + z + 1: #E = 116.
        let gf = Gf2::new(&IrreduciblePoly {
            degree: 7,
            low_terms: vec![0, 1],
        });
        assert_eq!(count_binary(&gf, 0, 1), 116);
        // y² = x³ + x + 13 over GF(59): #E = 67 (the 8-bit bench curve).
        assert_eq!(count_prime(59, 1, 13), 67);
    }
}
