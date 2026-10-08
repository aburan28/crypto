//! Disclosed n17 encoding/matrix diagnostic. No search or DLP solve.
use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, matrix_f4_f2_counted_with, permute_poly, FieldStructure,
    RowCriterion,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::polynomial_reuse::{
    build_decomposition_system_reusing, DecompositionTemplate,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use serde_json::{json, Value};

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
    let template = DecompositionTemplate::build(&fb.subspace_basis, &kc.curve.b, 3, &st).unwrap();
    let mut controls = Vec::new();
    for item in input["controls"].as_array().unwrap() {
        let x = item["target"][0].as_u64().unwrap();
        let positions: Vec<_> = (0..17).filter(|k| x & (1 << k) != 0).collect();
        let target = F2mElement::from_bit_positions(&positions, 17);
        let direct =
            build_decomposition_system(&fb.subspace_basis, &target, &kc.curve.b, 3, &st).unwrap();
        let cached =
            build_decomposition_system_reusing(&fb.subspace_basis, &target, &kc.curve.b, 3, &st)
                .unwrap();
        let repeated =
            build_decomposition_system_reusing(&fb.subspace_basis, &target, &kc.curve.b, 3, &st)
                .unwrap();
        let instantiated = template.instantiate(&target);
        let permutation = direct.interleaved_order(17);
        let interleaved: Vec<_> = cached
            .equations
            .iter()
            .map(|p| permute_poly(p, &permutation))
            .collect();
        let mut matrices = Vec::new();
        for (layout, polys) in [
            ("original", &cached.equations),
            ("interleaved", &interleaved),
        ] {
            for (engine, criterion) in [("f4", RowCriterion::None), ("f5", RowCriterion::F5)] {
                let result = matrix_f4_f2_counted_with(polys, direct.n_vars, 3, criterion);
                matrices.push(match result {
                    Some((reduced, words)) => json!({"layout":layout,"engine":engine,
                        "status":"REDUCED","rows":rows(&reduced),"elimination_and_criterion_word_xors":words}),
                    None => json!({"layout":layout,"engine":engine,"status":"OVERSIZE","rows":null,
                        "elimination_and_criterion_word_xors":null}),
                });
            }
        }
        controls.push(json!({"trial":item["trial"],"target":item["target"],
            "n_vars":direct.n_vars,"ell":direct.ell,"m":direct.m,
            "direct":rows(&direct.equations),"template":rows(&instantiated.equations),
            "reused":rows(&cached.equations),"repeated":rows(&repeated.equations),
            "permutation":permutation,"interleaved":rows(&interleaved),"matrices":matrices}));
    }
    println!(
        "{}",
        json!({"schema_version":1,"scope":"disclosed encoding and full-readback root matrices",
        "basis":fb.subspace_basis.iter().map(|e| e.raw_bits()[0]).collect::<Vec<_>>(),
        "controls":controls})
    );
}
