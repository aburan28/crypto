//! Guided disclosed specialization paths, not solver search or DLP recovery.
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, diagnostic_production_decisive, diagnostic_production_substitute,
    permute_poly, FieldStructure, RowCriterion,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use serde_json::{json, Value};
use std::collections::BTreeMap;

fn rows(polys: &[F2BoolPoly]) -> Vec<Vec<u64>> {
    polys
        .iter()
        .map(|p| p.terms.iter().map(|t| t.mask).collect())
        .collect()
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 2);
    let input: Value = serde_json::from_slice(&std::fs::read(&args[1]).unwrap()).unwrap();
    assert_eq!(input["fixture"]["degree"], 17);
    assert_eq!(input["fixture"]["curve_a"], 1);
    assert_eq!(input["summands"], 3);
    assert_eq!(input["matrix_degree"], 3);
    let kc = KoblitzCurve::new(1, 17).unwrap();
    let fb = build_standard_subspace_factor_base(&kc, 6).unwrap();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let depths = [0, 6, 12, 18, 24, 30, 35];
    let mut controls = Vec::new();
    for item in input["controls"].as_array().unwrap() {
        let x = item["target"][0].as_u64().unwrap();
        let positions: Vec<_> = (0..17).filter(|k| x & (1 << k) != 0).collect();
        let target = F2mElement::from_bit_positions(&positions, 17);
        let direct =
            build_decomposition_system(&fb.subspace_basis, &target, &kc.curve.b, 3, &st).unwrap();
        assert_eq!(direct.n_vars, 35);
        let permutation = direct.interleaved_order(17);
        let interleaved: Vec<_> = direct
            .equations
            .iter()
            .map(|p| permute_poly(p, &permutation))
            .collect();
        let mut paths = Vec::new();
        let mut nodes = Vec::new();
        let mut node_index = BTreeMap::new();
        for (layout, root) in [
            ("original", &direct.equations),
            ("interleaved", &interleaved),
        ] {
            for model in item["models"].as_array().unwrap() {
                let assignment = model[if layout == "original" {
                    "assignment"
                } else {
                    "renamed_assignment"
                }]
                .as_str()
                .unwrap()
                .parse::<u64>()
                .unwrap();
                let mut system = root.clone();
                let mut steps = Vec::new();
                let mut node_ids = Vec::new();
                for depth in 0..=35 {
                    steps.push(rows(&system));
                    if depths.contains(&depth) {
                        let key = serde_json::to_string(&rows(&system)).unwrap();
                        let id = if let Some(id) = node_index.get(&key) {
                            *id
                        } else {
                            let id = nodes.len();
                            let mut engines = Vec::new();
                            for (engine, criterion) in
                                [("f4", RowCriterion::None), ("f5", RowCriterion::F5)]
                            {
                                let result =
                                    diagnostic_production_decisive(&system, 35, 3, criterion);
                                engines.push(match result {
                                    Some(reduced) => json!({"engine":engine,"status":"REDUCED","rows":rows(&reduced)}),
                                    None => json!({"engine":engine,"status":"OVERSIZE","rows":null}),
                                });
                            }
                            nodes.push(json!({"system":rows(&system),"engines":engines}));
                            node_index.insert(key, id);
                            id
                        };
                        node_ids.push(id);
                    }
                    if depth < 35 {
                        system = system
                            .iter()
                            .map(|p| {
                                diagnostic_production_substitute(
                                    p,
                                    depth,
                                    assignment & (1 << depth) != 0,
                                )
                            })
                            .collect();
                    }
                }
                paths.push(json!({"layout":layout,"order":model["order"],"steps":steps,"node_ids":node_ids}));
            }
        }
        controls.push(
            json!({"trial":item["trial"],"n_vars":35,"permutation":permutation,
                            "paths":paths,"nodes":nodes}),
        );
    }
    println!(
        "{}",
        json!({"schema_version":1,"depths":depths,"controls":controls})
    );
}
