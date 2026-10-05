//! Frozen small-base candidate for research/f6_n83_xfilter_20261005.
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::f6_wide_geometry::F6SignedPairIndex;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use serde_json::json;

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 curve");
    let target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&target));
    assert_eq!(kc.mul(&target, &kc.subgroup_order), BinaryPoint::Infinity);
    for ell in [8u32, 10] {
        let parent = build_standard_subspace_factor_base(&kc, ell).unwrap();
        let base = cofactor_project_factor_base(&kc, &parent).unwrap();
        let b = base.points.len();
        let pairs = b * (b + 1) / 2;
        let warmup = F6SignedPairIndex::new(&kc.curve, &base.points, pairs).unwrap();
        let expected = warmup.solve4_xonly(&target);
        let signed_sums = warmup.signed_sum_count();
        drop(warmup);
        let mut build_ns = Vec::new();
        let mut query_ns = Vec::new();
        for _ in 0..3 {
            let start = Instant::now();
            let index = F6SignedPairIndex::new(&kc.curve, &base.points, pairs).unwrap();
            build_ns.push(start.elapsed().as_nanos());
            assert_eq!(index.pair_count(), pairs);
            assert_eq!(index.signed_sum_count(), signed_sums);
            let start = Instant::now();
            let answer = index.solve4_xonly(&target);
            query_ns.push(start.elapsed().as_nanos());
            assert_eq!(std::hint::black_box(answer.is_some()), expected.is_some());
            if let Some(indices) = answer {
                let replay = indices.iter().fold(BinaryPoint::Infinity, |sum, &i| {
                    kc.add(&sum, &base.points[i])
                });
                assert_eq!(replay, target);
            }
        }
        println!(
            "{}",
            json!({
                "curve_id": kc.label(),
                "ell": ell,
                "usable_points": b,
                "pairs": pairs,
                "signed_sums": signed_sums,
                "outcome": if expected.is_some() { "witness" } else { "no_witness" },
                "build_ns": build_ns,
                "query_ns": query_ns,
                "correctness": "PASS",
                "claim_scope": "exploratory_stage_only"
            })
        );
    }
}
