//! Bounded checks on the named curve's exact parameters. Public fixtures only;
//! this diagnostic does not collect a factor base or recover a logarithm.

use num_bigint::BigUint;
use num_traits::One;
use serde_json::{json, Value};

use crypto_lib::binary_ecc::curve::{point_add, scalar_mul, BinaryPoint};
use crypto_lib::cryptanalysis::binary_semaev::binary_semaev_s3;
use crypto_lib::cryptanalysis::curve_catalog::{CatalogCurve, CurveObject};
use crypto_lib::cryptanalysis::ec_index_calculus::semaev_s3;
use crypto_lib::ecc::point::Point;

fn hex(n: &BigUint) -> String {
    n.to_str_radix(16)
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
                            "s3_zero": semaev_s3(x1, x2, x3, &a, &c.b_fe()).is_zero(),
                            "wrong_target_rejected": p.add(&q, &a) != wrong,
                            "wrong_s3_rejected": !semaev_s3(x1, x2, bad, &a, &c.b_fe()).is_zero(),
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
    let verified = parameters_verified
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
        "full_parameter_pipeline_available": false,
        "scope": "Named-curve parameter, public point-sum and S3 checks only; no relation collection, independent-rank test, linear algebra or DLP recovery.",
    }))
}
