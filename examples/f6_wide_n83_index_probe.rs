//! Registered pair-index probe: research/f6_wide_n83_20261004/INDEX_PROTOCOL.md.
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::f6_wide_geometry::F6WidePairIndex;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use serde_json::json;

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned gate curve");
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
        let warmup = F6WidePairIndex::new(&kc.curve, &base.points, pairs).unwrap();
        let expected = warmup.solve4(&target);
        drop(warmup);
        let mut build_ns = Vec::new();
        let mut query_ns = Vec::new();
        for _ in 0..3 {
            let start = Instant::now();
            let index = F6WidePairIndex::new(&kc.curve, &base.points, pairs).unwrap();
            build_ns.push(start.elapsed().as_nanos());
            assert_eq!(index.pair_count(), pairs);
            let start = Instant::now();
            let answer = index.solve4(&target);
            query_ns.push(start.elapsed().as_nanos());
            assert_eq!(std::hint::black_box(answer), expected);
        }
        println!(
            "{}",
            json!({
                "curve_id": kc.label(),
                "ell": ell,
                "usable_points": b,
                "pairs": pairs,
                "outcome": if expected.is_some() { "witness" } else { "no_witness" },
                "build_ns": build_ns,
                "query_ns": query_ns,
                "correctness": "PASS",
                "claim_scope": "exploratory_stage_only"
            })
        );
    }
}
