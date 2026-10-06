use super::checker::*;
use num_bigint::BigUint;
use serde_json::{json, Value};

fn prime(p: u32, a: u32, b: u32) -> Value {
    json!({"v":"1","form":"y^2=x^3+a*x+b","p":p.to_string(),"field":format!("fp-{p}"),"a":a.to_string(),"b":b.to_string()})
}
fn binary(q: u32, a: u32, b: u32) -> Value {
    let m = 31 - q.leading_zeros();
    let hash = digest(format!("f2m-modulus:0x{q:x}").as_bytes());
    json!({"v":"1","form":"y^2+xy=x^3+a*x^2+b","modulus":format!("0x{q:x}"),"field":format!("f2m-{m}-{}",&hash[..8]),"a":a.to_string(),"b":b.to_string()})
}
fn checked(v: &Value) -> Model {
    match model(v) {
        Ok(m) => m,
        Err(_) => panic!("expected valid model"),
    }
}
fn check_affine_points(m: &Model, c: &Certificate, q: u32) {
    let k = &m.field;
    let h = decoded(&c.h, k).unwrap();
    let f = decoded(&c.f, k).unwrap();
    let x = decoded(&c.x, k).unwrap();
    let a = decoded(&c.y_v, k).unwrap();
    let b = decoded(&c.y_0, k).unwrap();
    for u in (0..q).map(BigUint::from) {
        for v in (0..q).map(BigUint::from) {
            if k.add(&k.mul(&v, &v), &k.mul(&k.eval(&h, &u), &v)) != k.eval(&f, &u) {
                continue;
            }
            let x = k.eval(&x, &u);
            let y = k.add(&k.mul(&k.eval(&a, &u), &v), &k.eval(&b, &u));
            let x2 = k.mul(&x, &x);
            let mut lhs = k.mul(&y, &y);
            let mut rhs = k.add(&k.mul(&x2, &x), &m.b);
            match k {
                Field::Prime(_) => rhs = k.add(&rhs, &k.mul(&m.a, &x)),
                Field::Binary(_) => {
                    lhs = k.add(&lhs, &k.mul(&x, &y));
                    rhs = k.add(&rhs, &k.mul(&m.a, &x2));
                }
            }
            assert_eq!(lhs, rhs);
        }
    }
}

#[test]
fn all_small_prime_models_and_cover_points() {
    for p in [5u32, 7, 11] {
        for a in 0..p {
            for b in 0..p {
                let v = prime(p, a, b);
                if (4 * a * a * a + 27 * b * b).is_multiple_of(p) {
                    assert!(model(&v).is_err());
                    continue;
                }
                let m = checked(&v);
                let c = construct(&m);
                verify(&m, &c).unwrap();
                if b == 0 {
                    assert_ne!(c.x[0], "0x0");
                }
                check_affine_points(&m, &c, p);
            }
        }
    }
}
#[test]
fn all_f8_models_and_cover_points() {
    for a in 0..8 {
        for b in 1..8 {
            let m = checked(&binary(0xb, a, b));
            let c = construct(&m);
            verify(&m, &c).unwrap();
            check_affine_points(&m, &c, 8);
        }
    }
    // Degree-one polynomial-basis field is also supported.
    let m = checked(&binary(3, 0, 1));
    verify(&m, &construct(&m)).unwrap();
}
#[test]
fn rejects_invalid_fields_singular_models_and_noncanonical_coefficients() {
    for v in [
        prime(5, 0, 0),
        prime(15, 1, 1),
        binary(0x15, 0, 1),
        binary(0xb, 1, 0),
        prime(5, 5, 1),
        binary(0xb, 8, 1),
    ] {
        assert!(matches!(model(&v), Err(ModelError::Invalid(_))));
    }
    let mut v = prime(7, 1, 1);
    v["field"] = json!("fp-11");
    assert!(matches!(model(&v), Err(ModelError::Invalid(_))));
    assert!(matches!(
        model(&prime(3, 1, 1)),
        Err(ModelError::Unsupported(_))
    ));
    // ICV1's extension part: listed coefficients over GF(p^k), unsupported
    // rather than invalid.
    let ext = json!({"v":"1","form":"y^2=x^3+a*x+b","p":"5","k":"2","modulus":["2","0"],
        "field":"fpk-5-2-90ac2fda","a":["1","0"],"b":["1","0"]});
    assert!(matches!(model(&ext), Err(ModelError::Unsupported(_))));
}
#[test]
fn rejects_tampered_certificates() {
    for v in [prime(11, 1, 1), binary(0xb, 1, 3)] {
        let m = checked(&v);
        let c = construct(&m);
        for field in ["genus", "degree"] {
            let mut bad = serde_json::to_value(&c).unwrap();
            bad[field] = json!(99);
            assert!(verify(&m, &serde_json::from_value(bad).unwrap()).is_err());
        }
        for field in ["f", "x", "y_0", "y_v", "h"] {
            let mut bad = serde_json::to_value(&c).unwrap();
            bad[field] = if bad[field].as_array().unwrap().is_empty() {
                json!(["0x1"])
            } else {
                json!([])
            };
            assert!(verify(&m, &serde_json::from_value(bad).unwrap()).is_err());
        }
        let mut bad = c.clone();
        bad.f.push("0x0".into());
        assert!(verify(&m, &bad).is_err());
        let mut bad = c.clone();
        bad.construction = "unknown".into();
        assert!(verify(&m, &bad).is_err());
    }
}
fn row(v: Value) -> Value {
    let raw = v.to_string();
    let hash = digest(raw.as_bytes());
    json!({"slug":format!("test-{}",&hash[..8]),"icv1":format!("ICV1:test:{}",&hash[..12]),"model_json":raw})
}
#[test]
fn unknown_is_not_no_cover_and_identities_are_bound() {
    let mut unknown = prime(7, 1, 1);
    unknown["form"] = json!("unknown");
    let mut input = json!({"schema_version":1,"curves":[row(unknown),row(prime(5,0,0))]});
    let report = catalog(input.to_string().as_bytes()).unwrap();
    assert_eq!(
        report["summary"],
        json!({"verified":0,"unsupported":1,"invalid_input":1})
    );
    for c in report["curves"].as_array().unwrap() {
        assert!(c["exists"].is_null());
    }
    let duplicate = input["curves"][0].clone();
    input["curves"].as_array_mut().unwrap().push(duplicate);
    assert!(catalog(input.to_string().as_bytes()).is_err());
    input["curves"].as_array_mut().unwrap().pop();
    input["curves"][0]["slug"] = json!("bad-hash");
    assert!(catalog(input.to_string().as_bytes()).is_err());
}
#[test]
fn entire_catalog_replays_deterministically() {
    let input = include_bytes!("../../../docs/curves/registry.json");
    let first = catalog(input).unwrap();
    assert_eq!(first, catalog(input).unwrap());
    let registry: Value = serde_json::from_slice(input).unwrap();
    // Every curve but the extension fields' (ICV1's extension part, B5a),
    // for which no construction here exists, has a verified cover.
    let rows = registry["curves"].as_array().unwrap();
    let extension = rows.iter().filter(|r| r["family"] == "extension").count();
    assert_eq!(
        first["summary"]["verified"].as_u64().unwrap() as usize,
        rows.len() - extension
    );
    assert_eq!(
        first["summary"]["unsupported"].as_u64().unwrap() as usize,
        extension
    );
    for (row, finding) in rows.iter().zip(first["curves"].as_array().unwrap()) {
        if row["family"] == "extension" {
            assert_eq!(finding["status"], "unsupported");
            assert!(finding["exists"].is_null());
            continue;
        }
        let m = checked(&serde_json::from_str(row["model_json"].as_str().unwrap()).unwrap());
        let c: Certificate = serde_json::from_value(finding["certificate"].clone()).unwrap();
        verify(&m, &c).unwrap();
        assert!(finding["dlp_advantage"].is_null());
        assert!(finding["minimal_genus"].is_null());
    }
}
