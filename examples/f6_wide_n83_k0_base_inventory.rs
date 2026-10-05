//! Structural gate inventory: research/f6_wide_n83_20261004/K0_BASE_PROTOCOL.md.
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use num_bigint::BigUint;
use serde_json::json;

fn choose4_with_replacement(b: usize) -> BigUint {
    let b = BigUint::from(b);
    (&b * (&b + 1u32) * (&b + 2u32) * (&b + 3u32)) / 24u32
}

fn main() {
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 gate failed validation");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    for ell in [8u32, 10, 12] {
        let parent = build_standard_subspace_factor_base(&kc, ell)
            .unwrap_or_else(|err| panic!("ell={ell} parent: {err}"));
        let base = cofactor_project_factor_base(&kc, &parent)
            .unwrap_or_else(|err| panic!("ell={ell} projection: {err}"));
        let b = base.points.len();
        println!(
            "{}",
            json!({
                "curve_id": kc.label(),
                "ell": ell,
                "parent_points": parent.points.len(),
                "usable_points": b,
                "signed_columns": base.signed_orbits.len(),
                "unordered4": choose4_with_replacement(b).to_string(),
                "subgroup_order": kc.subgroup_order.to_string(),
                "status": "PASS"
            })
        );
    }
}
