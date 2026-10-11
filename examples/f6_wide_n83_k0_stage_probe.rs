//! Frozen K0 stage diagnostic: research/f6_wide_n83_20261004/K0_STAGE_PROTOCOL.md.
use std::time::Instant;

use crypto_lib::binary_ecc::curve::point_add;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::f6_wide_geometry::{batch_add_fixed, F6WidePairIndex};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use serde_json::json;

fn scalar_pairs(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Vec<BinaryPoint> {
    let mut sums = Vec::with_capacity(points.len() * (points.len() + 1) / 2);
    for (i, fixed) in points.iter().enumerate() {
        for other in &points[i..] {
            sums.push(point_add(&kc.curve, fixed, other));
        }
    }
    sums
}

fn batch_pairs(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Vec<BinaryPoint> {
    let mut sums = Vec::with_capacity(points.len() * (points.len() + 1) / 2);
    for (i, fixed) in points.iter().enumerate() {
        sums.extend(batch_add_fixed(&kc.curve, fixed, &points[i..]));
    }
    sums
}

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned gate curve");
    let parent = build_standard_subspace_factor_base(&kc, 8).expect("ell8 base");
    let base = cofactor_project_factor_base(&kc, &parent).expect("usable projection");
    let points = &base.points[..64];
    let target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&target));
    assert_eq!(kc.mul(&target, &kc.subgroup_order), BinaryPoint::Infinity);

    let expected = scalar_pairs(&kc, points);
    assert_eq!(expected.len(), 2080);
    assert_eq!(batch_pairs(&kc, points), expected);
    let _ = std::hint::black_box(scalar_pairs(&kc, points));
    let _ = std::hint::black_box(batch_pairs(&kc, points));
    let mut scalar_pair_ns = Vec::new();
    let mut batch_pair_ns = Vec::new();
    for round in 0..3 {
        for method in if round % 2 == 0 { [0, 1] } else { [1, 0] } {
            let start = Instant::now();
            let result = if method == 0 {
                scalar_pairs(&kc, points)
            } else {
                batch_pairs(&kc, points)
            };
            let elapsed = start.elapsed().as_nanos();
            assert_eq!(std::hint::black_box(result), expected);
            if method == 0 {
                scalar_pair_ns.push(elapsed);
            } else {
                batch_pair_ns.push(elapsed);
            }
        }
    }

    let index = F6WidePairIndex::new(&kc.curve, points, 2080).expect("capped index");
    let expected_witness = index.solve4_scalar_reference(&target);
    assert_eq!(index.solve4(&target), expected_witness);
    let _ = std::hint::black_box(index.solve4_scalar_reference(&target));
    let _ = std::hint::black_box(index.solve4(&target));
    let mut scalar_query_ns = Vec::new();
    let mut batch_query_ns = Vec::new();
    for round in 0..3 {
        for method in if round % 2 == 0 { [0, 1] } else { [1, 0] } {
            let start = Instant::now();
            let result = if method == 0 {
                index.solve4_scalar_reference(&target)
            } else {
                index.solve4(&target)
            };
            let elapsed = start.elapsed().as_nanos();
            assert_eq!(std::hint::black_box(result), expected_witness);
            if method == 0 {
                scalar_query_ns.push(elapsed);
            } else {
                batch_query_ns.push(elapsed);
            }
        }
    }
    println!(
        "{}",
        json!({
            "curve_id": kc.label(),
            "target_x": "0x355fb5df7a905f16921eb",
            "target_y": "0x5900a390f42d290f1bbe",
            "ell": 8,
            "usable_base_points": base.points.len(),
            "subset_points": points.len(),
            "pairs": index.pair_count(),
            "outcome": if expected_witness.is_some() { "witness" } else { "no_witness" },
            "scalar_pair_ns": scalar_pair_ns,
            "batch_pair_ns": batch_pair_ns,
            "scalar_query_ns": scalar_query_ns,
            "batch_query_ns": batch_query_ns,
            "correctness": "PASS",
            "claim_scope": "exploratory_stage_only"
        })
    );
}
