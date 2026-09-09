//! Finite checks of the two-torsion quotient coordinate u=x+1/x.
//! No unknown scalar recovery. Counts rational fibers and verifies endomorphism identities.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{points_with_x, KoblitzCurve};
use num_bigint::BigUint;
use serde_json::json;
use std::collections::HashMap;

fn trace(x: &F2mElement, kc: &KoblitzCurve) -> bool {
    let mut t = F2mElement::zero(kc.n);
    let mut q = x.clone();
    for _ in 0..kc.n {
        t = t.add(&q);
        q = q.square(&kc.curve.irreducible);
    }
    assert!(t.is_zero() || t == F2mElement::one(kc.n));
    !t.is_zero()
}
fn quotient(p: &BinaryPoint, kc: &KoblitzCurve) -> BinaryPoint {
    let mut q = kc.add(p, p);
    for _ in 1..kc.n {
        q = kc.frobenius(&q);
    }
    q
}
fn audit(a: u8, n: u32) -> serde_json::Value {
    let kc = KoblitzCurve::new(a, n).unwrap();
    let mut fibers = HashMap::new();
    let mut points = 0;
    for x in 0..(1u32 << n) {
        let x = F2mElement::from_biguint(&BigUint::from(x), n);
        for p in points_with_x(&kc.curve, &x) {
            points += 1;
            let q = quotient(&p, &kc);
            if x.is_zero() {
                assert_eq!(q, BinaryPoint::Infinity);
                continue;
            }
            let u = x.add(&x.flt_inverse(&kc.curve.irreducible).unwrap());
            let BinaryPoint::Affine { x: qx, .. } = &q else {
                panic!("non-torsion point maps to infinity")
            };
            assert_eq!(qx, &u);
            *fibers.entry(u.to_biguint()).or_insert(0usize) += 1;
            assert_eq!(
                quotient(&kc.add(&p, kc.generator()), &kc),
                kc.add(&q, &quotient(kc.generator(), &kc))
            );
            if kc.cofactor == BigUint::from(2u8) {
                assert_eq!(kc.mul(&q, &kc.subgroup_order), BinaryPoint::Infinity);
            }
        }
    }
    let mut admissible = 0;
    for u in 1..(1u32 << n) {
        let u = F2mElement::from_biguint(&BigUint::from(u), n);
        let expected = trace(&u, &kc) == (a == 1 && n % 2 == 1)
            && !trace(&u.flt_inverse(&kc.curve.irreducible).unwrap(), &kc);
        let count = fibers.get(&u.to_biguint()).copied().unwrap_or(0);
        assert_eq!(count, if expected { 4 } else { 0 });
        admissible += usize::from(expected);
    }
    let zero_fiber = fibers.get(&BigUint::from(0u8)).copied().unwrap_or(0);
    assert_eq!(zero_fiber, if a == 0 { 2 } else { 0 }); // odd n
    json!({"kind":"quotient_geometry_audit","a":a,"n":n,"affine_points":points,
        "nonzero_u_values":(1u32<<n)-1,"admissible_nonzero_u":admissible,"generic_fiber_size":4,
        "u_zero_fiber_size":zero_fiber,"endomorphism_coordinate_checks":points-1,
        "cofactor":kc.cofactor.to_string(),"scope":"rational fibers, trace predicates and group identities"})
}
fn main() {
    for (a, n) in [(0, 7), (1, 7), (0, 9), (1, 9)] {
        println!("{}", audit(a, n));
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quotient_trace_predicate_matches_complete_fibers() {
        audit(0, 7);
        audit(1, 7);
    }
}
