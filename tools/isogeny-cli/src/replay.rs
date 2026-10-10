//! Independent exact certificate replay, compiled into the installed tool.
use super::{
    curve,
    field::{Fe, Field},
    kernel, map_identity,
    poly::{self, Poly},
};
use crate::sha256_hex;
use curve::Model;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::{fs, path::Path};
fn json_number(n: &BigUint) -> Value {
    serde_json::from_str(&n.to_string()).unwrap()
}
fn integer(v: &Value) -> BigUint {
    let s = v.as_str().expect("decimal or hex integer string");
    if let Some(h) = s.strip_prefix("0x") {
        BigUint::parse_bytes(h.as_bytes(), 16).unwrap()
    } else {
        BigUint::parse_bytes(s.as_bytes(), 10).unwrap()
    }
}
fn coeffs(f: &Field, v: &Value) -> Poly {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| f.from_big(&integer(x)))
        .collect()
}
fn eval(f: &Field, n: &Poly, d: &Poly, p: Option<(Fe, Fe)>) -> Option<(Fe, Fe)> {
    let (x, y) = p?;
    let nv = poly::eval(f, n, &x);
    let dv = poly::eval(f, d, &x);
    if f.is_zero(&dv) {
        return None;
    }
    let top = f.sub(
        &f.mul(&poly::eval(f, &poly::derivative(f, n), &x), &dv),
        &f.mul(&nv, &poly::eval(f, &poly::derivative(f, d), &x)),
    );
    Some((
        f.div(&nv, &dv).unwrap(),
        f.mul(&y, &f.div(&top, &f.sqr(&dv)).unwrap()),
    ))
}
pub fn verify(out: &Path, write_certificate: bool) -> Value {
    let registry: Value = serde_json::from_str(crate::REGISTRY).unwrap();
    let summary: Value =
        serde_json::from_slice(&fs::read(out.join("search.json")).unwrap()).unwrap();
    let output_suffix = "";
    let mut coverage = None;
    let mut replayed = vec![];
    let mut curves = vec![];
    for attempt in summary["attempts"].as_array().unwrap() {
        let preset = attempt["curve"].as_str().unwrap();
        let ell = attempt["ell"].as_u64().unwrap();
        assert_eq!(attempt["receipt"]["stdout"], format!("ell-{ell}.json"));
        let path = out
            .join(preset)
            .join(attempt["receipt"]["stdout"].as_str().unwrap());
        let bytes = fs::read(&path).unwrap();
        assert_eq!(
            sha256_hex(&bytes),
            attempt["receipt"]["stdout_sha256"].as_str().unwrap()
        );
        let stderr = fs::read(out.join(preset).join(format!("ell-{ell}.stderr.txt"))).unwrap();
        assert_eq!(
            sha256_hex(&stderr),
            attempt["receipt"]["stderr_sha256"].as_str().unwrap()
        );
        if attempt["receipt"]["status"] != "PASS" {
            continue;
        }
        assert_eq!(attempt["receipt"]["exit_code"], 0);
        let doc: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(doc["status"], "PASS");
        assert_eq!(doc["curve"]["preset"], preset);
        let p = integer(&doc["curve"]["p"]);
        let aa = integer(&doc["curve"]["a"]);
        let b = integer(&doc["curve"]["b"]);
        let f = Field::new(&p).unwrap();
        let source = Model {
            a: f.from_big(&aa),
            b: f.from_big(&b),
        };
        let source_slug = doc["result"]["source_icv1"]["slug"].as_str().unwrap();
        let source_entry = registry["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["slug"] == source_slug)
            .unwrap();
        assert_eq!(integer(&source_entry["params"]["p"]), p);
        assert_eq!(integer(&source_entry["params"]["a"]), aa);
        assert_eq!(integer(&source_entry["params"]["b"]), b);
        let group_order = integer(&source_entry["order"]);
        let standard =
            crate::catalog::representation(source_entry).expect("no complete source subgroup");
        let n = integer(&standard["curve"]["subgroup_order"]);
        let cofactor = crate::catalog::integer(&standard["curve"]["cofactor"]);
        let cofactor_big = BigUint::parse_bytes(cofactor.to_string().as_bytes(), 10).unwrap();
        assert_eq!(
            num_bigint::BigUint::parse_bytes(cofactor.to_string().as_bytes(), 10).unwrap() * &n,
            group_order
        );
        let g = (
            f.from_big(&integer(&standard["curve"]["generator"][0])),
            f.from_big(&integer(&standard["curve"]["generator"][1])),
        );
        assert!(source.on_curve(&f, &g.0, &g.1));
        assert!(curve::scalar_mul(&f, &source, &g.0, &g.1, &n).is_none());
        let isos = doc["result"]["isogenies"].as_array().unwrap();
        let partial = doc["scope"] == "one_frobenius_eigenline";
        let expected = if partial { "PARTIAL" } else { "COMPLETE" };
        if partial {
            assert_eq!(doc["degree_coverage"], "PARTIAL");
            assert_eq!(doc["expected_total_eigenlines"], 2);
            assert_eq!(doc["constructed_eigenlines"], 1);
            assert_eq!(isos.len(), 1);
        } else {
            assert_eq!(doc["scope"], "all_frobenius_eigenlines");
            assert_eq!(doc["degree_coverage"], "COMPLETE");
            assert_eq!(isos.len(), 2);
        }
        if let Some(previous) = coverage {
            assert_eq!(previous, expected);
        }
        coverage = Some(expected);
        for (index, iso) in isos.iter().enumerate() {
            assert_eq!(iso["verified"], true);
            let h = coeffs(&f, &iso["kernel"]);
            let num = coeffs(&f, &iso["x_map"]["num"]);
            let den = coeffs(&f, &iso["x_map"]["den"]);
            let cod = kernel::verify_kernel(&f, &source, ell, &h)
                .expect("independent kernel certificate failed");
            assert_eq!(f.to_big(&cod.a), integer(&iso["codomain"]["a"]));
            assert_eq!(f.to_big(&cod.b), integer(&iso["codomain"]["b"]));
            assert_eq!(
                poly::mul(&f, &h, &h),
                den,
                "denominator must equal kernel squared"
            );
            assert_eq!(num.len(), ell as usize + 1);
            assert_eq!(num.last(), Some(&f.one()));
            assert_eq!(poly::degree(&poly::gcd(&f, &num, &den)), Some(0));
            assert!(
                map_identity::check(&f, &source, &cod, &num, &den),
                "rational map fails exact curve-equation substitution"
            );
            let gp = eval(&f, &num, &den, Some(g))
                .expect("prime subgroup generator cannot lie in the small-degree kernel");
            assert!(cod.on_curve(&f, &gp.0, &gp.1));
            assert!(curve::scalar_mul(&f, &cod, &gp.0, &gp.1, &n).is_none());
            let mut point_checks = 0;
            for seed_x in [101, 1009, 10007, 65537] {
                let q = curve::least_point(&f, &source, seed_x);
                let iq = eval(&f, &num, &den, Some(q)).unwrap();
                assert!(cod.on_curve(&f, &iq.0, &iq.1));
                for k in [2, 3, 17, 65537, ell] {
                    let k = BigUint::from(k);
                    let left = eval(
                        &f,
                        &num,
                        &den,
                        curve::scalar_mul(&f, &source, &q.0, &q.1, &k),
                    );
                    let right = curve::scalar_mul(&f, &cod, &iq.0, &iq.1, &k);
                    assert_eq!(left, right, "public scalar transport failed");
                    point_checks += 1;
                }
            }
            assert!(eval(&f, &num, &den, None).is_none());
            let target_id = &iso["codomain"]["icv1"];
            let target_a = f.to_big(&cod.a);
            let target_b = f.to_big(&cod.b);
            let field_record = json!({"characteristic":json_number(&p),"degree":1,"representation":"prime","element_encoding":"hex integer modulo the characteristic"});
            let curve_record = json!({"model":"short Weierstrass","coefficients":[0,0,0,json_number(&target_a),json_number(&target_b)],"subgroup_order":n.to_string(),"cofactor":json_number(&cofactor_big),
                "generator":[format!("0x{:x}",f.to_big(&gp.0)),format!("0x{:x}",f.to_big(&gp.1))],"target_group":"prime-order subgroup"});
            let digest = sha256_hex(
                serde_json::to_string(&json!({"field":field_record,"curve":curve_record}))
                    .unwrap()
                    .as_bytes(),
            );
            let ec1 = format!("EC1P{}Cfph{}", p.bits(), &digest[..12]);
            let uid = format!("urn:ec-record:1:sha256:{digest}");
            let row = json!({"source_icv1":source_slug,"source_ec1":standard["ec1"],"source_curve_uid":standard["curve_uid"],
                "ell":ell,"index":index,"target_icv1":target_id["slug"],"target_ec1":ec1,"target_curve_uid":uid,
                "kernel_degree":h.len()-1,"map_numerator_coefficients":num.len(),"map_denominator_coefficients":den.len(),
                "certificate_sha256":sha256_hex(&bytes),"kernel_check":"PASS","codomain_check":"PASS","subgroup_check":"PASS","exact_rational_map_check":"PASS",
                "public_scalar_transport_checks":point_checks,"status":"PASS"});
            replayed.push(row);
            curves.push(json!({"name":target_id["slug"],"icv1":target_id["icv1"],"p":p.to_string(),"a":target_a.to_string(),"b":target_b.to_string(),
                "group_order":group_order.to_string(),"subgroup_order":n.to_string(),"cofactor":json_number(&cofactor_big),"j":iso["codomain"]["j"],"generator":[format!("0x{:x}",f.to_big(&gp.0)),format!("0x{:x}",f.to_big(&gp.1))],
                "source":source_slug,"ell":ell,"index":index,"target_ec1":ec1,"target_curve_uid":uid,"model_json":target_id["model_json"]}));
        }
        eprintln!(
            "replayed {preset} ell={ell}: {} certified maps; coverage {}",
            isos.len(),
            coverage.unwrap()
        );
    }
    assert!(!replayed.is_empty(), "no completed construction to certify");
    let receipt = json!({"schema":"large-degree-isogeny-replay/v1","method":"existing walker kernel verifier, separate field and polynomial implementation",
        "scope":"same-host independent implementation","degree_coverage":coverage.unwrap(),"expected_total_eigenlines":2,"records":replayed,"status":"PASS"});
    if write_certificate {
        fs::write(
            out.join(format!("replay{output_suffix}.json")),
            format!("{}\n", serde_json::to_string_pretty(&receipt).unwrap()),
        )
        .unwrap();
        fs::write(
            out.join(format!("curves{output_suffix}.json")),
            format!(
                "{}\n",
                serde_json::to_string_pretty(
                    &json!({"schema":"large-degree-isogeny-curves/v1","curves":curves})
                )
                .unwrap()
            ),
        )
        .unwrap();
    }
    receipt
}
