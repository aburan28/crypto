//! Exact subgroup-coordinate support census for the frozen different-j bases.

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ext_curve::{ExtCurve, ExtPoint, Fp2};
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps};
use serde_json::{json, Value};

type Point = ExtPoint<[u64; 2]>;
type Curve = ExtCurve<Fp2>;

fn parse_point(v: &Value) -> Point {
    if v["infinity"].as_bool().unwrap() {
        return Point::infinity([0, 0]);
    }
    let coords = |name: &str| -> [u64; 2] {
        let pair = v[name].as_array().unwrap();
        [pair[0].as_u64().unwrap(), pair[1].as_u64().unwrap()]
    };
    Point::affine(coords("x"), coords("y"))
}

fn point_list(v: &Value) -> Vec<Point> {
    v.as_array().unwrap().iter().map(parse_point).collect()
}

fn base_support(
    curve: &Curve,
    index: &BTreeMap<u64, u64>,
    base: &[Point],
) -> (Vec<u64>, BTreeSet<u64>) {
    assert_eq!(base.len(), 8);
    let mut base_scalars = Vec::new();
    let mut signed_orbits = BTreeSet::new();
    for &point in base {
        assert!(curve.is_on_curve(point));
        let scalar = index[&curve.key(&point)];
        assert!(scalar > 0 && scalar < 61);
        assert!(signed_orbits.insert(scalar.min(61 - scalar)));
        base_scalars.push(scalar);
    }
    let signed: Vec<u64> = base_scalars
        .iter()
        .flat_map(|&scalar| [scalar, 61 - scalar])
        .collect();
    let mut support = BTreeSet::new();
    for i in 0..signed.len() {
        for j in i..signed.len() {
            let sum = (signed[i] + signed[j]) % 61;
            if sum != 0 {
                support.insert(sum);
            }
        }
    }
    (base_scalars, support)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 2, "usage: different_j_support_census PUBLIC_V2");
    let raw = fs::read(&args[1]).unwrap();
    let public: Value = serde_json::from_slice(&raw).unwrap();
    assert_eq!(public["schema"], "different_j_public/v2");
    let field = Fp2::new(67).unwrap();
    let source = Curve {
        f: field,
        a: [46, 0],
        b: [37, 0],
    };
    let destination = Curve {
        f: field,
        a: [53, 0],
        b: [25, 0],
    };
    let generator = parse_point(&public["generator"]);
    let image_generator = parse_point(&public["mapped_generator"]);
    let mut ops = GroupOps::default();
    let mut source_index = BTreeMap::new();
    let mut destination_index = BTreeMap::new();
    for scalar in 0..61 {
        let source_point = source.mul(&mut ops, generator, scalar);
        let image_point = destination.mul(&mut ops, image_generator, scalar);
        assert!(source_index
            .insert(source.key(&source_point), scalar)
            .is_none());
        assert!(destination_index
            .insert(destination.key(&image_point), scalar)
            .is_none());
    }
    assert_eq!(source_index.len(), 61);
    assert_eq!(destination_index.len(), 61);
    let bases = [
        ("source", &source, &source_index, "base_projected"),
        (
            "isogeny-image",
            &destination,
            &destination_index,
            "mapped_base_projected",
        ),
        (
            "destination-direct",
            &destination,
            &destination_index,
            "direct_base_projected",
        ),
        (
            "full-image",
            &destination,
            &destination_index,
            "full_base_projected",
        ),
    ];
    let mut rows = BTreeMap::new();
    let mut supports = BTreeMap::new();
    for (name, curve, index, key) in bases {
        let (scalars, support) = base_support(curve, index, &point_list(&public[key]));
        rows.insert(
            name,
            json!({"base_scalars":scalars,"support_size":support.len(),"support":support}),
        );
        supports.insert(name, support);
    }
    assert_eq!(supports["source"], supports["isogeny-image"]);
    let full = &supports["full-image"];
    let projected = &supports["isogeny-image"];
    let direct = &supports["destination-direct"];
    println!(
        "{}",
        json!({
            "schema":"different_j_support_census/v1",
            "public_sha256":sha256_hex(&raw),
            "r":61,
            "arms":rows,
            "full_vs_projected":{
                "overlap":full.intersection(projected).count(),
                "full_only":full.difference(projected).copied().collect::<Vec<_>>(),
                "projected_only":projected.difference(full).copied().collect::<Vec<_>>()
            },
            "full_vs_destination_direct":{
                "overlap":full.intersection(direct).count(),
                "full_only":full.difference(direct).copied().collect::<Vec<_>>(),
                "direct_only":direct.difference(full).copied().collect::<Vec<_>>()
            },
            "group_ops":ops
        })
    );
}
