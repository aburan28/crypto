//! Frozen stage diagnostic: research/f6_wide_n83_20261004/PROTOCOL.md.
use std::time::Instant;

use crypto_lib::binary_ecc::curve::{point_add, scalar_mul};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::f6_wide_geometry::{batch_add_fixed, F6WidePairIndex};
use num_bigint::BigUint;
use serde_json::json;

fn curve() -> BinaryCurve {
    let m = 83;
    let parse = |s: &str| BigUint::parse_bytes(s.as_bytes(), 16).unwrap();
    BinaryCurve {
        m,
        irreducible: IrreduciblePoly {
            degree: m,
            low_terms: vec![0, 1, 2, 45],
        },
        a: F2mElement::one(m),
        b: F2mElement::one(m),
        generator: BinaryPoint::Affine {
            x: F2mElement::from_biguint(&parse("68a212cfe19a809fe0598"), m),
            y: F2mElement::from_biguint(&parse("244a245ea0b17d8cc8297"), m),
        },
        order: BigUint::from(8_569_786_107_849_059u64),
        cofactor: BigUint::from(1_128_547_018u64),
    }
}

fn scalar_pairs(curve: &BinaryCurve, points: &[BinaryPoint]) -> Vec<BinaryPoint> {
    let mut sums = Vec::with_capacity(points.len() * (points.len() + 1) / 2);
    for (i, fixed) in points.iter().enumerate() {
        for other in &points[i..] {
            sums.push(point_add(curve, fixed, other));
        }
    }
    sums
}

fn batch_pairs(curve: &BinaryCurve, points: &[BinaryPoint]) -> Vec<BinaryPoint> {
    let mut sums = Vec::with_capacity(points.len() * (points.len() + 1) / 2);
    for (i, fixed) in points.iter().enumerate() {
        sums.extend(batch_add_fixed(curve, fixed, &points[i..]));
    }
    sums
}

fn main() {
    let curve = curve();
    assert!(curve.is_on_curve(&curve.generator));
    assert_eq!(
        scalar_mul(&curve, &curve.generator, &curve.order),
        BinaryPoint::Infinity
    );
    let points: Vec<_> = (1u32..=64)
        .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
        .collect();
    for size in [16usize, 32, 64] {
        let base = &points[..size];
        let expected = scalar_pairs(&curve, base);
        assert_eq!(batch_pairs(&curve, base), expected);
        let _ = std::hint::black_box(scalar_pairs(&curve, base));
        let _ = std::hint::black_box(batch_pairs(&curve, base));
        let mut scalar_ns = Vec::new();
        let mut batch_ns = Vec::new();
        for round in 0..3 {
            let order = if round % 2 == 0 { [0, 1] } else { [1, 0] };
            for method in order {
                let start = Instant::now();
                let result = if method == 0 {
                    scalar_pairs(&curve, base)
                } else {
                    batch_pairs(&curve, base)
                };
                let elapsed = start.elapsed().as_nanos();
                assert_eq!(std::hint::black_box(result), expected);
                if method == 0 {
                    scalar_ns.push(elapsed);
                } else {
                    batch_ns.push(elapsed);
                }
            }
        }
        let index = F6WidePairIndex::new(&curve, base, expected.len()).unwrap();
        let planted = scalar_mul(&curve, &curve.generator, &BigUint::from(2 * size as u32));
        assert!(index.solve4(&planted).is_some());
        let missing = scalar_mul(
            &curve,
            &curve.generator,
            &BigUint::from(4 * size as u32 + 1),
        );
        assert!(index.solve4(&missing).is_none());
        assert!(index.solve4_scalar_reference(&missing).is_none());
        let _ = std::hint::black_box(index.solve4(&missing));
        let _ = std::hint::black_box(index.solve4_scalar_reference(&missing));
        let mut scalar_query_ns = Vec::new();
        let mut batch_query_ns = Vec::new();
        for round in 0..3 {
            let order = if round % 2 == 0 { [0, 1] } else { [1, 0] };
            for method in order {
                let start = Instant::now();
                let result = if method == 0 {
                    index.solve4_scalar_reference(&missing)
                } else {
                    index.solve4(&missing)
                };
                let elapsed = start.elapsed().as_nanos();
                assert!(std::hint::black_box(result).is_none());
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
                "n": 83,
                "subgroup_order": curve.order.to_string(),
                "size": size,
                "pairs": index.pair_count(),
                "reference_ns": scalar_ns,
                "batch_ns": batch_ns,
                "scalar_query_ns": scalar_query_ns,
                "batch_query_ns": batch_query_ns,
                "correctness": "PASS",
                "claim_scope": "exploratory_stage_only"
            })
        );
    }
}
