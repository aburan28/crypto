//! Planted full-base n=83 K0 witness control; not a natural-yield sample.
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::f6_wide_geometry::F6SignedPairIndex;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use serde_json::json;

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 curve");
    let parent = build_standard_subspace_factor_base(&kc, 12).unwrap();
    let base = cofactor_project_factor_base(&kc, &parent).unwrap();
    assert_eq!(base.points.len(), 4_054);
    let planted_indices = [0usize, 2, 4, 6];
    let target = planted_indices
        .into_iter()
        .fold(BinaryPoint::Infinity, |sum, i| {
            kc.add(&sum, &base.points[i])
        });
    assert_ne!(target, BinaryPoint::Infinity);
    assert!(kc.curve.is_on_curve(&target));
    let pairs = base.points.len() * (base.points.len() + 1) / 2;
    let index = F6SignedPairIndex::new(&kc.curve, &base.points, pairs).unwrap();
    let reference = index.solve4_xonly(&target).expect("reference witness");
    let packed = index.solve4_xonly_pmull83(&target).expect("packed witness");
    for indices in [reference, packed] {
        let replay = indices.iter().fold(BinaryPoint::Infinity, |sum, &i| {
            kc.add(&sum, &base.points[i])
        });
        assert_eq!(replay, target);
    }
    println!(
        "{}",
        json!({
            "curve_id": kc.label(),
            "usable_points": base.points.len(),
            "pairs": pairs,
            "planted_indices": planted_indices,
            "reference_witness": reference,
            "candidate_witness": packed,
            "correctness": "PASS",
            "claim_scope": "planted_control_not_natural_yield"
        })
    );
}
