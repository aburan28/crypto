//! Structural screen: research/f6_n83_mixed_base_20261004/ORBIT_PROTOCOL.md.
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use num_bigint::BigUint;
use serde_json::json;
use std::collections::HashSet;

fn choose_with_replacement(b: usize, k: u32) -> BigUint {
    let mut count = BigUint::from(1u8);
    for i in 0..k {
        count *= BigUint::from(b + i as usize);
        count /= BigUint::from(i + 1);
    }
    count
}

fn point_bytes(point: &BinaryPoint) -> Vec<u8> {
    match point {
        BinaryPoint::Infinity => vec![0],
        BinaryPoint::Affine { x, y } => {
            let mut bytes = Vec::with_capacity(33);
            bytes.push(1);
            for word in x.raw_bits().iter().chain(y.raw_bits()) {
                bytes.extend_from_slice(&word.to_le_bytes());
            }
            bytes
        }
    }
}

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 curve failed validation");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    assert_eq!(
        kc.frobenius(kc.generator()),
        kc.mul(kc.generator(), &kc.lambda)
    );
    let geometric = build_standard_subspace_factor_base(&kc, 12).expect("geometric seed base");
    let projected = cofactor_project_factor_base(&kc, &geometric).expect("projected seed base");
    assert_eq!(projected.points.len(), 4054);
    assert_eq!(projected.signed_orbits.len(), 2027);
    let mut closure = HashSet::new();
    let mut canonical_orbits = HashSet::new();
    for seed in &projected.points {
        let mut current = seed.clone();
        let mut canonical: Option<Vec<u8>> = None;
        for _ in 0..kc.extension_degree() {
            let positive = point_bytes(&current);
            let negative = point_bytes(&point_neg(&current));
            for candidate in [positive, negative] {
                if canonical
                    .as_ref()
                    .is_none_or(|best| candidate.as_slice() < best.as_slice())
                {
                    canonical = Some(candidate.clone());
                }
                closure.insert(candidate);
            }
            current = kc.frobenius(&current);
        }
        assert_eq!(&current, seed, "Frobenius orbit did not close at 83");
        canonical_orbits.insert(canonical.expect("nonempty orbit"));
    }
    assert_eq!(
        closure.len(),
        canonical_orbits.len() * 2 * kc.extension_degree() as usize
    );
    assert!(canonical_orbits.len() <= projected.signed_orbits.len());
    let mut ordered: Vec<_> = closure.iter().collect();
    ordered.sort_unstable();
    let mut hasher = blake3::Hasher::new();
    for point in ordered {
        hasher.update(point);
    }
    let b = closure.len();
    println!(
        "{}",
        json!({
            "curve_id":kc.label(),
            "seed_ell":12,
            "seed_usable_points":projected.points.len(),
            "seed_signed_columns":projected.signed_orbits.len(),
            "frobenius_order":kc.extension_degree(),
            "closure_usable_points":b,
            "closure_signed_frobenius_columns":canonical_orbits.len(),
            "closure_set_blake3":hasher.finalize().to_hex().to_string(),
            "unordered4_ceiling_numerator":choose_with_replacement(b,4).to_string(),
            "unordered5_ceiling_numerator":choose_with_replacement(b,5).to_string(),
            "subgroup_order":kc.subgroup_order.to_string(),
            "status":"PASS"
        })
    );
}
