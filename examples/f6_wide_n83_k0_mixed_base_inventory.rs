//! Structural screen: research/f6_n83_mixed_base_20261004/PROTOCOL.md.
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
    let subgroup_order = &kc.subgroup_order;
    let mut terminal: Option<HashSet<Vec<u8>>> = None;
    let mut terminal_b = 0usize;
    for ell in [12u32, 14, 15, 16] {
        let result = build_standard_subspace_factor_base(&kc, ell).and_then(|parent| {
            cofactor_project_factor_base(&kc, &parent).map(|base| (parent, base))
        });
        let (parent, base) = match result {
            Ok(pair) => pair,
            Err(error) => {
                println!("{}", json!({"ell":ell,"status":"ERROR","error":error}));
                continue;
            }
        };
        let points: HashSet<_> = base.points.iter().map(point_bytes).collect();
        let mut ordered: Vec<_> = points.iter().collect();
        ordered.sort_unstable();
        let mut hasher = blake3::Hasher::new();
        for point in ordered {
            hasher.update(point);
        }
        let b = base.points.len();
        let six = choose_with_replacement(b, 6);
        let subset = terminal.as_ref().map(|set| set.is_subset(&points));
        let mixed = if terminal_b > 0 {
            Some(choose_with_replacement(b, 4) * choose_with_replacement(terminal_b, 2))
        } else {
            None
        };
        println!(
            "{}",
            json!({
                "curve_id":kc.label(),
                "ell":ell,
                "parent_points":parent.points.len(),
                "usable_points":b,
                "signed_columns":base.signed_orbits.len(),
                "projected_set_blake3":hasher.finalize().to_hex().to_string(),
                "terminal_subset":subset,
                "unordered6_ceiling_numerator":six.to_string(),
                "mixed4plus2_ceiling_numerator":mixed.map(|value|value.to_string()),
                "subgroup_order":subgroup_order.to_string(),
                "status":"PASS"
            })
        );
        if ell == 12 {
            terminal_b = b;
            terminal = Some(points);
        }
    }
}
