//! Read-only inventory for the n9 F6-IC cold-solve control.
//! This is not part of the timed worker. It exposes actual usable points
//! after cofactor projection so candidate IDs never use a nominal bound.
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_factor_base_search::FactorBaseSpec;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    cofactor_project_factor_base, projected_signed_orbit_count, KoblitzCurve,
};
use serde_json::json;

fn main() {
    let curve = KoblitzCurve::new(1, 9).expect("n9 a1 curve");
    let geometric = FactorBaseSpec::StandardSubspace { dimension: 4 }
        .materialize(&curve)
        .expect("standard subspace base");
    let usable =
        cofactor_project_factor_base(&curve, &geometric).expect("nonempty subgroup-usable base");
    let identity_images = geometric
        .points
        .iter()
        .filter(|point| curve.mul(point, &curve.cofactor) == BinaryPoint::Infinity)
        .count();
    let encode = |points: &[crypto_lib::binary_ecc::BinaryPoint]| {
        let mut out: Vec<_> = points
            .iter()
            .map(|point| {
                let BinaryPoint::Affine { x, y } = point else {
                    panic!("factor base must not contain identity")
                };
                [x.to_biguint().to_string(), y.to_biguint().to_string()]
            })
            .collect();
        out.sort_by(|a, b| {
            let ax = a[0].parse::<u64>().expect("small field");
            let ay = a[1].parse::<u64>().expect("small field");
            let bx = b[0].parse::<u64>().expect("small field");
            let by = b[1].parse::<u64>().expect("small field");
            (ax, ay).cmp(&(bx, by))
        });
        out
    };
    println!(
        "{}",
        json!({
            "schema_version":1,
            "curve_a":1,
            "degree":9,
            "recipe":{"kind":"standard_subspace","dimension":4},
            "geometric_points":encode(&geometric.points),
            "usable_points":encode(&usable.points),
            "identity_images":identity_images,
            "duplicate_nonidentity_images":geometric.points.len()-identity_images-usable.points.len(),
            "effective_columns":projected_signed_orbit_count(&curve,&geometric)
        })
    );
}
