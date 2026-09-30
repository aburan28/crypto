//! Bounded checks on the named curve's exact parameters. Public fixtures only;
//! this diagnostic does not collect a factor base or recover a logarithm.

use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};

use crypto_lib::binary_ecc::curve::{
    point_add, point_double, point_neg, scalar_mul, BinaryCurve, BinaryPoint,
};
use crypto_lib::binary_ecc::f2m::F2mElement;
use crypto_lib::cryptanalysis::binary_semaev::{
    binary_semaev_s3, solve_artin_schreier, solve_quadratic_f2m,
};
use crypto_lib::cryptanalysis::curve_catalog::{CatalogCurve, CurveObject};
use crypto_lib::cryptanalysis::ec_index_calculus::semaev_s3;
use crypto_lib::ecc::point::Point;

fn hex(n: &BigUint) -> String {
    n.to_str_radix(16)
}

/// Exercise point lifting and exceptional group operations on public inputs.
fn binary_diagnostics(c: &BinaryCurve) -> Result<Value, String> {
    let BinaryPoint::Affine { x, y } = &c.generator else {
        return Err("binary diagnostics require an affine generator".into());
    };
    let irr = &c.irreducible;
    let zero = F2mElement::zero(c.m);
    let one = F2mElement::one(c.m);
    let rhs = x
        .square(irr)
        .mul(x, irr)
        .add(&c.a.mul(&x.square(irr), irr))
        .add(&c.b);
    let roots = solve_quadratic_f2m(&one, x, &rhs, c.m, irr);
    let mut encoded_roots: Vec<String> = roots.iter().map(|r| hex(&r.to_biguint())).collect();
    encoded_roots.sort_unstable();
    let square_roots = solve_quadratic_f2m(&one, &zero, &c.b, c.m, irr);
    let torsion_verified = square_roots.len() == 1 && {
        let torsion = BinaryPoint::Affine {
            x: zero.clone(),
            y: square_roots[0].clone(),
        };
        c.is_on_curve(&torsion)
            && point_double(c, &torsion) == BinaryPoint::Infinity
            && point_add(c, &torsion, &torsion) == BinaryPoint::Infinity
    };
    let artin_input = y.square(irr).add(y);
    let artin_verified = solve_artin_schreier(&artin_input, c.m, irr)
        .is_some_and(|r| r.square(irr).add(&r) == artin_input);
    let g = &c.generator;
    Ok(json!({
        "generator_y_roots": encoded_roots,
        "checks": {
            "generator_lift_two_roots": roots.len() == 2 && roots[0] != roots[1],
            "generator_lift_matches_signed_points": roots.contains(y) && roots.contains(&y.add(x)),
            "generator_lift_roots_on_curve": roots.iter().all(|r| c.is_on_curve(&BinaryPoint::Affine { x: x.clone(), y: r.clone() })),
            "linear_root_verified": solve_quadratic_f2m(&zero, &one, y, c.m, irr) == vec![y.clone()],
            "nonzero_constant_has_no_roots": solve_quadratic_f2m(&zero, &zero, &one, c.m, irr).is_empty(),
            "pure_square_and_two_torsion_verified": torsion_verified,
            "artin_schreier_round_trip": artin_verified,
            "artin_schreier_zero": solve_artin_schreier(&zero, c.m, irr) == Some(zero),
            "trace_one_equation_classified": solve_artin_schreier(&one, c.m, irr).is_none() == (c.m % 2 == 1),
            "left_identity": point_add(c, &BinaryPoint::Infinity, g) == *g,
            "right_identity": point_add(c, g, &BinaryPoint::Infinity) == *g,
            "inverse_sum": point_add(c, g, &point_neg(g)) == BinaryPoint::Infinity,
            "double_matches_sum": point_double(c, g) == point_add(c, g, g),
            "zero_scalar": scalar_mul(c, g, &BigUint::zero()) == BinaryPoint::Infinity,
            "order_minus_one_matches_inverse": scalar_mul(c, g, &(&c.order - BigUint::one())) == point_neg(g),
        },
    }))
}

/// Export every parameter used by the diagnostic, without narrowing to u64.
pub fn exact_parameters(curve: &CatalogCurve) -> Value {
    match &curve.object {
        CurveObject::Prime(c) => json!({
            "model": "short-weierstrass", "encoding": "hex-integer",
            "field": {"characteristic": hex(&c.p), "degree": 1},
            "a": hex(&c.a), "b": hex(&c.b),
            "generator": {"x": hex(&c.gx), "y": hex(&c.gy)},
            "subgroup_order": hex(&c.n), "cofactor": hex(&BigUint::from(c.h)),
        }),
        CurveObject::Binary(c) => {
            let mut terms = c.irreducible.low_terms.clone();
            terms.push(c.irreducible.degree);
            terms.sort_unstable();
            let generator = match &c.generator {
                BinaryPoint::Affine { x, y } => {
                    json!({"x": hex(&x.to_biguint()), "y": hex(&y.to_biguint())})
                }
                BinaryPoint::Infinity => Value::Null,
            };
            json!({
                "model": "binary-weierstrass", "encoding": "hex-polynomial-bits",
                "field": {"characteristic": "2", "degree": c.m,
                    "basis": "polynomial", "modulus_exponents": terms},
                "a": hex(&c.a.to_biguint()), "b": hex(&c.b.to_biguint()),
                "generator": generator,
                "subgroup_order": hex(&c.order), "cofactor": hex(&c.cofactor),
            })
        }
        CurveObject::Char3(_) => Value::Null,
    }
}

/// Check known point sums and S3, including public scalars wider than a word.
pub fn validate(curve: &CatalogCurve) -> Result<Value, String> {
    if matches!(&curve.object, CurveObject::Char3(_)) {
        return Err("full-parameter IC diagnostics currently cover prime and binary curves".into());
    }
    let checks = curve.verify();
    let parameters_verified = checks.iter().all(|c| c.passed);
    let checks_json: Vec<Value> = checks
        .iter()
        .map(|c| json!({"name": c.name, "passed": c.passed, "detail": c.detail}))
        .collect();
    let order = curve.subgroup_order();
    if order <= BigUint::from(7u32) {
        return Err("full-parameter fixtures require a subgroup order greater than seven".into());
    }
    let scalars = [
        (BigUint::one(), BigUint::from(2u32)),
        (&order - BigUint::from(2u32), &order - BigUint::from(3u32)),
    ];
    let mut fixtures = Vec::new();
    // Invalid parameters must never reach unchecked affine inversions.
    if parameters_verified {
        for (u, v) in scalars {
            let w = (&u + &v) % &order;
            let wrong_w = (&w + BigUint::one()) % &order;
            let mut fixture = match &curve.object {
                CurveObject::Prime(c) => {
                    let a = c.a_fe();
                    let b = c.fe(c.b.clone());
                    let p = c.generator().scalar_mul(&u, &a);
                    let q = c.generator().scalar_mul(&v, &a);
                    let r = c.generator().scalar_mul(&w, &a);
                    let wrong = c.generator().scalar_mul(&wrong_w, &a);
                    match (&p, &q, &r, &wrong) {
                        (
                            Point::Affine { x: x1, .. },
                            Point::Affine { x: x2, .. },
                            Point::Affine { x: x3, y: y3 },
                            Point::Affine { x: bad, .. },
                        ) => json!({
                            "target": {"x": hex(&x3.value), "y": hex(&y3.value)},
                            "point_sum_verified": p.add(&q, &a) == r,
                            "s3_zero": semaev_s3(x1, x2, x3, &a, &b).is_zero(),
                            "wrong_target_rejected": p.add(&q, &a) != wrong,
                            "wrong_s3_rejected": !semaev_s3(x1, x2, bad, &a, &b).is_zero(),
                        }),
                        _ => return Err("unexpected infinity in prime IC fixture".into()),
                    }
                }
                CurveObject::Binary(c) => {
                    let p = scalar_mul(c, &c.generator, &u);
                    let q = scalar_mul(c, &c.generator, &v);
                    let r = scalar_mul(c, &c.generator, &w);
                    let wrong = scalar_mul(c, &c.generator, &wrong_w);
                    match (&p, &q, &r, &wrong) {
                        (
                            BinaryPoint::Affine { x: x1, .. },
                            BinaryPoint::Affine { x: x2, .. },
                            BinaryPoint::Affine { x: x3, y: y3 },
                            BinaryPoint::Affine { x: bad, .. },
                        ) => json!({
                            "target": {"x": hex(&x3.to_biguint()), "y": hex(&y3.to_biguint())},
                            "point_sum_verified": point_add(c, &p, &q) == r,
                            "s3_zero": binary_semaev_s3(x1, x2, x3, &c.b, &c.irreducible).is_zero(),
                            "wrong_target_rejected": point_add(c, &p, &q) != wrong,
                            "wrong_s3_rejected": !binary_semaev_s3(x1, x2, bad, &c.b, &c.irreducible).is_zero(),
                        }),
                        _ => return Err("unexpected infinity in binary IC fixture".into()),
                    }
                }
                CurveObject::Char3(_) => unreachable!("guarded above"),
            };
            fixture["u"] = json!(hex(&u));
            fixture["v"] = json!(hex(&v));
            fixture["w"] = json!(hex(&w));
            fixtures.push(fixture);
        }
    }
    let binary_diagnostics = if parameters_verified {
        match &curve.object {
            CurveObject::Binary(c) => binary_diagnostics(c)?,
            _ => Value::Null,
        }
    } else {
        Value::Null
    };
    let binary_verified = binary_diagnostics.is_null()
        || binary_diagnostics["checks"]
            .as_object()
            .is_some_and(|checks| checks.values().all(|v| v == true));
    let verified = parameters_verified
        && binary_verified
        && fixtures.len() == 2
        && fixtures.iter().all(|f| {
            [
                "point_sum_verified",
                "s3_zero",
                "wrong_target_rejected",
                "wrong_s3_rejected",
            ]
            .iter()
            .all(|k| f[*k] == true)
        });
    Ok(json!({
        "schema_version": 1, "operation": "validate",
        "status": if verified { "checks_passed" } else { "checks_failed" },
        "verified": verified, "curve": curve.name, "family": curve.family.tag(),
        "field_bits": curve.field().bits, "order_bits": curve.order_bits(),
        "is_named_curve": true, "diagnostic_only": true,
        "exact_parameters": exact_parameters(curve), "checks": checks_json,
        "fixtures": fixtures,
        "binary_diagnostics": binary_diagnostics,
        "full_parameter_pipeline_available": false,
        "scope": "Named-curve parameter, public point-sum and S3 checks only; no relation collection, independent-rank test, linear algebra or DLP recovery.",
    }))
}
