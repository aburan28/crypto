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
/// A model over `GF(p)[t]/(t^k + low)`, coefficients as digit lists.
fn ext(p: u32, low: &[u32], a: &[u32], b: &[u32]) -> Value {
    let list = |d: &[u32]| d.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    let joined = list(low).join(",");
    let hash = digest(format!("fpk-modulus:{p}:{joined}").as_bytes());
    json!({"v":"1","form":"y^2=x^3+a*x+b","p":p.to_string(),"k":low.len().to_string(),
        "modulus":list(low),"field":format!("fpk-{p}-{}-{}",low.len(),&hash[..8]),
        "a":list(a),"b":list(b)})
}
/// The base-`p` digits of `n`, `k` of them.
fn digits(p: u32, k: usize, mut n: u32) -> Vec<u32> {
    (0..k)
        .map(|_| {
            let d = n % p;
            n /= p;
            d
        })
        .collect()
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
        let (hu, fu) = (k.eval(&h, &u), k.eval(&f, &u));
        for v in (0..q).map(BigUint::from) {
            if k.add(&k.mul(&v, &v), &k.mul(&hu, &v)) != fu {
                continue;
            }
            let x = k.eval(&x, &u);
            let y = k.add(&k.mul(&k.eval(&a, &u), &v), &k.eval(&b, &u));
            let x2 = k.mul(&x, &x);
            let mut lhs = k.mul(&y, &y);
            let mut rhs = k.add(&k.mul(&x2, &x), &m.b);
            match k {
                Field::Prime(_) | Field::Extension { .. } => rhs = k.add(&rhs, &k.mul(&m.a, &x)),
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
fn all_gf25_models_and_cover_points() {
    // GF(25) = GF(5)[t]/(t^2 + 2): -2 = 3 is not a square mod 5.
    let mut verified = 0;
    for a in 0..25 {
        for b in 0..25 {
            let v = ext(5, &[2, 0], &digits(5, 2, a), &digits(5, 2, b));
            let Ok(m) = model(&v) else {
                // Singular: 4a^3 + 27b^2 = 0 in GF(25).
                let k = Field::Extension {
                    p: BigUint::from(5u32),
                    low: vec![BigUint::from(2u32), BigUint::from(0u32)],
                };
                let (a, b) = (BigUint::from(a), BigUint::from(b));
                let disc = k.add(
                    &k.mul(&k.constant(4), &k.mul(&k.mul(&a, &a), &a)),
                    &k.mul(&k.constant(27), &k.mul(&b, &b)),
                );
                assert_eq!(disc, BigUint::from(0u32));
                continue;
            };
            let c = construct(&m);
            verify(&m, &c).unwrap();
            check_affine_points(&m, &c, 25);
            verified += 1;
        }
    }
    // 4a^3 = -27b^2 has 25 solutions (a, b) over GF(25): the singular models.
    assert_eq!(verified, 625 - 25);
}
#[test]
fn larger_extension_models_verify() {
    // GF(49) = GF(7)[t]/(t^2 + 1), every model; GF(125) = GF(5)[t]/(t^3 +
    // t + 1), a few models with every cover point.
    for a in 0..49 {
        for b in 0..49 {
            if let Ok(m) = model(&ext(7, &[1, 0], &digits(7, 2, a), &digits(7, 2, b))) {
                verify(&m, &construct(&m)).unwrap();
            }
        }
    }
    for (a, b) in [(0, 1), (1, 0), (5, 7), (124, 61)] {
        let m = checked(&ext(5, &[1, 1, 0], &digits(5, 3, a), &digits(5, 3, b)));
        let c = construct(&m);
        verify(&m, &c).unwrap();
        check_affine_points(&m, &c, 125);
    }
}
#[test]
fn rabin_counts_the_irreducible_polynomials() {
    // Monic irreducibles of degree k over GF(p): (1/k) sum_{d|k} mu(d) p^(k/d).
    for (p, k, expected) in [
        (5u32, 2usize, 10usize),
        (5, 3, 40),
        (5, 4, 150),
        (7, 2, 21),
        (2, 4, 3),
        (3, 3, 8),
    ] {
        let count = (0..p.pow(k as u32))
            .filter(|n| {
                let low: Vec<BigUint> = digits(p, k, *n).into_iter().map(BigUint::from).collect();
                irreducible_over_prime(&BigUint::from(p), &low)
            })
            .count();
        assert_eq!(count, expected, "p = {p}, k = {k}");
    }
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
    // ICV1's extension part: ICV1.md's own example verifies, ...
    let example = json!({"v":"1","form":"y^2=x^3+a*x+b","p":"5","k":"2","modulus":["2","0"],
        "field":"fpk-5-2-90ac2fda","a":["1","0"],"b":["1","0"]});
    assert_eq!(example, ext(5, &[2, 0], &[1, 0], &[1, 0]));
    let m = checked(&example);
    verify(&m, &construct(&m)).unwrap();
    // ... and a reducible modulus (t^2 - 1), a field label that names
    // another modulus, a coefficient not reduced mod p, a list of the wrong
    // length, a singular model and a degree below 2 are invalid input.
    let mut relabelled = ext(5, &[2, 0], &[1, 0], &[1, 0]);
    relabelled["field"] = ext(5, &[3, 0], &[1, 0], &[1, 0])["field"].clone();
    let mut short = ext(5, &[2, 0], &[1, 0], &[1, 0]);
    short["a"] = json!(["1"]);
    let mut degree_one = ext(5, &[2, 0], &[1, 0], &[1, 0]);
    degree_one["k"] = json!("1");
    for v in [
        ext(5, &[4, 0], &[1, 0], &[1, 0]),
        relabelled,
        ext(5, &[2, 0], &[5, 0], &[1, 0]),
        short,
        ext(5, &[2, 0], &[0, 0], &[0, 0]),
        degree_one,
    ] {
        assert!(matches!(model(&v), Err(ModelError::Invalid(_))));
    }
    // Characteristic 3 has no construction here.
    assert!(matches!(
        model(&ext(3, &[1, 0], &[1, 0], &[1, 0])),
        Err(ModelError::Unsupported(_))
    ));
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
    // Every curve, the extension fields' included, has a verified cover.
    let rows = registry["curves"].as_array().unwrap();
    assert!(rows.iter().any(|r| r["family"] == "extension"));
    assert_eq!(
        first["summary"]["verified"].as_u64().unwrap() as usize,
        rows.len()
    );
    for (row, finding) in rows.iter().zip(first["curves"].as_array().unwrap()) {
        let m = checked(&serde_json::from_str(row["model_json"].as_str().unwrap()).unwrap());
        let c: Certificate = serde_json::from_value(finding["certificate"].clone()).unwrap();
        verify(&m, &c).unwrap();
        assert!(finding["dlp_advantage"].is_null());
        assert!(finding["minimal_genus"].is_null());
    }
}

#[test]
fn all_small_model_changes_have_rational_inverses() {
    use crypto_lib::utils::mod_inverse;
    use num_traits::ToPrimitive;
    for p in [5u32, 7, 11] {
        let inv = |x: u32| {
            mod_inverse(&BigUint::from(x), &BigUint::from(p))
                .unwrap()
                .to_u32()
                .unwrap()
        };
        for form in [
            "B*y^2=x^3+A*x^2+x",
            "a*x^2+y^2=1+d*x^2*y^2",
            "x^2+y^2=c^2*(1+d*x^2*y^2)",
        ] {
            for a in 0..p {
                for b in 0..p {
                    let v = json!({"v":"1","field":format!("fp-{p}"),"p":p.to_string(),"form":form,"A":a.to_string(),"B":b.to_string(),"a":a.to_string(),"c":a.to_string(),"d":b.to_string()});
                    let singular = if form.starts_with("B*") {
                        b == 0 || a * a % p == 4 % p
                    } else if form.starts_with("a*") {
                        a == 0 || b == 0 || a == b
                    } else {
                        a == 0 || b == 0 || b * a * a * a * a % p == 1
                    };
                    if singular {
                        assert!(model(&v).is_err());
                        continue;
                    }
                    let m = checked(&v);
                    let n = m.normalization.as_ref().unwrap();
                    let get =
                        |key: &str| number(n[key].as_str().unwrap()).unwrap().to_u32().unwrap();
                    let mb = get("montgomery_B");
                    let shift = get("shift_A_over_3");
                    let scale = get("edwards_scale");
                    let sa = m.a.to_u32().unwrap();
                    let sb = m.b.to_u32().unwrap();
                    for x in 0..p {
                        for y in 0..p {
                            if y * y % p != (x * x * x + sa * x + sb) % p {
                                continue;
                            }
                            let u = (mb * x + p - shift) % p;
                            let w = mb * y % p;
                            let (tx, ty) = if form.starts_with("B*") {
                                (u, w)
                            } else {
                                if w == 0 || (u + 1) % p == 0 {
                                    continue;
                                }
                                (
                                    scale * u * inv(w) % p,
                                    scale * (u + p - 1) * inv((u + 1) % p) % p,
                                )
                            };
                            if form.starts_with("B*") {
                                assert_eq!(b * ty * ty % p, (tx * tx * tx + a * tx * tx + tx) % p);
                            } else if form.starts_with("a*") {
                                assert_eq!(
                                    (a * tx * tx + ty * ty) % p,
                                    (1 + b * tx * tx * ty * ty) % p
                                );
                            } else {
                                assert_eq!(
                                    (tx * tx + ty * ty) % p,
                                    a * a * (1 + b * tx * tx * ty * ty) % p
                                );
                            }
                            let (ru, rw) = if form.starts_with("B*") {
                                (tx, ty)
                            } else {
                                if tx == 0 || (scale + p - ty) % p == 0 {
                                    continue;
                                }
                                let ru = (scale + ty) * inv((scale + p - ty) % p) % p;
                                (ru, scale * ru * inv(tx) % p)
                            };
                            assert_eq!(((ru + shift) * inv(mb) % p, rw * inv(mb) % p), (x, y));
                        }
                    }
                    let c = construct(&m);
                    verify(&m, &c).unwrap();
                    let mut altered = c.clone();
                    altered.target_model_map.as_mut().unwrap()["montgomery_B"] = json!("0x0");
                    assert!(verify(&m, &altered).is_err());
                }
            }
        }
    }
}

#[test]
fn cover_graph_is_content_addressed_and_rejects_wrong_ec1_binding() {
    let v = prime(5, 0, 1);
    let raw = v.to_string();
    let mut rep = json!({"field":{"characteristic":5,"degree":1,"representation":"prime"},"curve":{"model":"short Weierstrass","coefficients":[0,0,0,0,1],"subgroup_order":"3","cofactor":2,"generator":["0x0","0x1"]}});
    fn identify(rep: &mut Value) {
        let h = digest(
            json!({"field":rep["field"],"curve":rep["curve"]})
                .to_string()
                .as_bytes(),
        );
        rep["ec1"] = json!(format!("EC1P3Ctesth{}", &h[..12]));
        rep["curve_uid"] = json!(format!("urn:ec-record:1:sha256:{h}"));
    }
    identify(&mut rep);
    let hash = digest(raw.as_bytes());
    let mut registry = json!({"schema_version":1,"curves":[{"slug":format!("test-{}",&hash[..8]),"icv1":format!("test:{}",&hash[..12]),"model_json":raw,"order":"6","representations":[rep]}]});
    let report = catalog(registry.to_string().as_bytes()).unwrap();
    let graph = super::links::graph(&registry, &report).unwrap();
    for (collection, key, prefix) in [
        ("covers", "cover_uid", "urn:hc-model:1:sha256:"),
        ("maps", "map_uid", "urn:curve-cover-map:1:sha256:"),
    ] {
        for record in graph[collection].as_object().unwrap().values() {
            assert_eq!(
                record[key],
                format!(
                    "{prefix}{}",
                    digest(record["record"].to_string().as_bytes())
                )
            );
        }
    }
    registry["curves"][0]["representations"][0]["curve"]["coefficients"][4] = json!(2);
    identify(&mut registry["curves"][0]["representations"][0]);
    assert!(super::links::graph(&registry, &report).is_err());
}
