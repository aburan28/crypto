use super::checker::{self, ModelError};
use crypto_lib::{cryptanalysis::curve_id, utils::mod_inverse};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
type R<T> = Result<T, String>;
fn num(v: &BigUint) -> Value {
    serde_json::from_str(&v.to_string()).expect("unsigned integer")
}
fn raw(v: &Value) -> R<BigUint> {
    checker::number(v.as_str().ok_or("missing parameter string")?)
}
fn parameter(v: &Value, k: &str) -> R<BigUint> {
    raw(&v["params"][k]["raw"])
}
fn prime_parameter(v: &Value, k: &str, p: &BigUint) -> R<BigUint> {
    let s = v["params"][k]["raw"]
        .as_str()
        .ok_or("missing coefficient")?;
    if let Some(s) = s.strip_prefix('-') {
        Ok((p - checker::number(s)? % p) % p)
    } else {
        Ok(checker::number(s)? % p)
    }
}
fn model_error(e: ModelError) -> String {
    match e {
        ModelError::Invalid(s) | ModelError::Unsupported(s) => s,
    }
}

fn import(v: &Value) -> R<Value> {
    let r = raw(&v["order"])?;
    let h = raw(&v["cofactor"])?;
    if r <= BigUint::one() || h.is_zero() {
        return Err("invalid subgroup/cofactor metadata".into());
    }
    let order = &r * &h;
    let mut representation = None;
    let (id, family, params, model) = match v["field"]["type"].as_str() {
        Some("Prime") => {
            let p = raw(&v["field"]["p"])?;
            if p <= BigUint::from(3u32) || p.bits() > 4096 {
                return Err("unsupported characteristic/field size".into());
            }
            let mut model = json!({"v":"1","field":format!("fp-{p}"),"p":p.to_string()});
            let keys: Vec<&str> = match v["form"].as_str() {
                Some("Weierstrass") => {
                    model["form"] = json!("y^2=x^3+a*x+b");
                    vec!["a", "b"]
                }
                Some("Montgomery") => {
                    model["form"] = json!("B*y^2=x^3+A*x^2+x");
                    vec!["a", "b"]
                }
                Some("TwistedEdwards") => {
                    model["form"] = json!("a*x^2+y^2=1+d*x^2*y^2");
                    vec!["a", "d"]
                }
                Some("Edwards") => {
                    model["form"] = json!("x^2+y^2=c^2*(1+d*x^2*y^2)");
                    vec!["c", "d"]
                }
                _ => return Err("unsupported prime model form".into()),
            };
            for key in &keys {
                let x = prime_parameter(v, key, &p)?;
                let k = if v["form"] == "Montgomery" {
                    key.to_uppercase()
                } else {
                    key.to_string()
                };
                model[k] = json!(x.to_string());
            }
            let m = checker::model(&model).map_err(model_error)?;
            let mut id =
                curve_id::prime(&p, &m.a, &m.b, &order).ok_or("order outside Hasse interval")?;
            if m.normalization.is_some() {
                let raw = model.to_string();
                let hash = checker::digest(raw.as_bytes());
                id.icv1 = format!("{}:{}", id.icv1.rsplit_once(':').unwrap().0, &hash[..12]);
                id.slug = format!("{}-{}", id.slug.rsplit_once('-').unwrap().0, &hash[..8]);
                id.model_json = raw;
            }
            let field = json!({"characteristic":num(&p),"degree":1,"representation":"prime","element_encoding":"hex integer modulo the characteristic"});
            if let (Ok(gx), Ok(gy)) = (
                raw(&v["generator"]["x"]["raw"]),
                raw(&v["generator"]["y"]["raw"]),
            ) {
                if gx >= p || gy >= p {
                    return Err("generator is not canonically encoded".into());
                }
                let (x, y) = if let Some(n) = &m.normalization {
                    let b = raw(&n["montgomery_B"])?;
                    let bi = mod_inverse(&b, &p).ok_or("invalid normalization")?;
                    let shift = raw(&n["shift_A_over_3"])?;
                    let (u, w) = if n["kind"] == "montgomery_from_short" {
                        (gx.clone(), gy.clone())
                    } else {
                        let c = raw(&n["edwards_scale"])?;
                        let u = (&c + &gy)
                            * mod_inverse(&((&c + &p - &gy) % &p), &p)
                                .ok_or("exceptional generator")?
                            % &p;
                        let w = &c * &u * mod_inverse(&gx, &p).ok_or("exceptional generator")? % &p;
                        (u, w)
                    };
                    ((u + shift) * &bi % &p, w * bi % &p)
                } else {
                    (gx.clone(), gy.clone())
                };
                if &y * &y % &p != (&x * &x * &x + &m.a * &x + &m.b) % &p {
                    return Err("generator fails the curve equation".into());
                }
                let coefficients = if v["form"] == "Weierstrass" {
                    json!([0, 0, 0, num(&m.a), num(&m.b)])
                } else {
                    let mut o = serde_json::Map::new();
                    for key in &keys {
                        o.insert(key.to_string(), num(&prime_parameter(v, key, &p)?));
                    }
                    Value::Object(o)
                };
                let curve = json!({"model":if v["form"]=="Weierstrass"{"short Weierstrass"}else{v["form"].as_str().unwrap()},"coefficients":coefficients,
                    "subgroup_order":r.to_string(),"cofactor":num(&h),"generator":[format!("0x{gx:x}"),format!("0x{gy:x}")],"target_group":"prime-order subgroup"});
                representation = Some(rep(field, curve, &format!("P{}", p.bits())));
            }
            let params =
                json!({"p":p.to_string(),"a":model["a"],"b":model["b"],"model_form":model["form"]});
            (id, "prime", params, model)
        }
        Some("Binary") => {
            if v["field"]["basis"] != "poly" {
                return Err(
                    "normal-basis coordinates require a sourced, verified basis conversion".into(),
                );
            }
            let degree = v["field"]["degree"]
                .as_u64()
                .ok_or("missing binary degree")?;
            if degree == 0 || degree > 4096 {
                return Err("invalid binary degree".into());
            }
            let mut modulus = BigUint::zero();
            let mut exponents = vec![];
            for term in v["field"]["poly"]
                .as_array()
                .ok_or("missing defining polynomial")?
            {
                let power = term["power"].as_u64().ok_or("missing polynomial power")?;
                if power > degree || raw(&term["coeff"])? != BigUint::one() {
                    return Err("invalid binary modulus term".into());
                }
                if modulus.bit(power) {
                    return Err("duplicate modulus term".into());
                }
                modulus.set_bit(power, true);
                exponents.push(power);
            }
            exponents.sort_by(|a, b| b.cmp(a));
            let a = parameter(v, "a")?;
            let b = parameter(v, "b")?;
            let koblitz = b.is_one() && (a.is_zero() || a.is_one());
            let id = curve_id::binary(
                degree as u32,
                &modulus,
                &a,
                &b,
                &order,
                if koblitz { Some(-7) } else { None },
            )
            .ok_or("invalid binary curve/order")?;
            let model: Value = serde_json::from_str(&id.model_json).map_err(|e| e.to_string())?;
            // Validate unreduced source coefficients too; do not let an identity constructor silently reduce them.
            if a.bits() > degree || b.bits() > degree {
                return Err("noncanonical binary coefficients".into());
            }
            let m = checker::model(&model).map_err(model_error)?;
            if let (Ok(x), Ok(y)) = (
                raw(&v["generator"]["x"]["raw"]),
                raw(&v["generator"]["y"]["raw"]),
            ) {
                if x.bits() > degree || y.bits() > degree {
                    return Err(
                        "compressed or noncanonical generator; full coordinates required".into(),
                    );
                }
                let k = &m.field;
                let x2 = k.mul(&x, &x);
                if k.add(&k.mul(&y, &y), &k.mul(&x, &y))
                    != k.add(&k.add(&k.mul(&x2, &x), &k.mul(&a, &x2)), &b)
                {
                    return Err("binary generator fails the curve equation".into());
                }
                let field = json!({"characteristic":2,"degree":degree,"representation":"polynomial","modulus_exponents":exponents,"element_encoding":"hex polynomial coefficient bitset, least significant bit is constant"});
                let curve = json!({"model":"binary Weierstrass","coefficients":[1,num(&a),0,0,num(&b)],"subgroup_order":r.to_string(),"cofactor":num(&h),"generator":[format!("0x{x:x}"),format!("0x{y:x}")],"target_group":"prime-order subgroup"});
                representation = Some(rep(field, curve, &format!("N{degree}")));
            }
            let params = if koblitz {
                json!({"n":degree,"a":if a.is_zero(){0}else{1},"modulus":format!("0x{modulus:x}")})
            } else {
                json!({"m":degree,"a":format!("0x{a:x}"),"b":format!("0x{b:x}"),"modulus":format!("0x{modulus:x}")})
            };
            (
                id,
                if koblitz { "koblitz" } else { "binary" },
                params,
                model,
            )
        }
        _ => {
            return Err(
                "extension/tower field adapter not implemented; exact source record retained"
                    .into(),
            )
        }
    };
    let m = checker::model(&model).map_err(model_error)?;
    let cert = checker::construct(&m);
    checker::verify(&m, &cert)?;
    let parts: Vec<_> = id.icv1.split(':').collect();
    let category = v["category"].as_str().unwrap_or("unknown");
    let adopted = !matches!(category, "other" | "bn" | "bls" | "mnt" | "nums");
    let source_id = v["source_id"].as_str().ok_or("missing source identifier")?;
    let mut names = vec![source_id.to_string()];
    if adopted {
        names.push(v["name"].as_str().ok_or("missing name")?.to_string());
    }
    let reps = representation.into_iter().collect::<Vec<_>>();
    Ok(
        json!({"slug":id.slug,"icv1":id.icv1,"model_json":id.model_json,"family":family,"params":params,
        "trace":serde_json::from_str::<Value>(parts[2]).map_err(|e|e.to_string())?,"order":parts[3],"j":parts[4],"end":parts[5],
        "aliases":names,"standard_names":if adopted{vec![v["name"].as_str().unwrap().to_string()]}else{vec![]},
        "sources":["docs/curves/standards/registry.json"],"representations":reps,
        "standards_provenance":[{"source_id":source_id,"display_name":v["name"],"category":category,"oid":v["oid"],"sources":v["sources"],
            "adoption_status":if adopted{"listed_by_cited_standard_or_standard_example"}else{"not_asserted"},
            "parameter_validation":"field and nonsingularity checked; supplied full generator checked on curve; group order remains source-declared"}]}),
    )
}

fn rep(field: Value, curve: Value, prefix: &str) -> Value {
    let hash = checker::digest(json!({"field":field,"curve":curve}).to_string().as_bytes());
    json!({"ec1":format!("EC1{prefix}Cstandardh{}",&hash[..12]),"curve_uid":format!("urn:ec-record:1:sha256:{hash}"),"field_sha256":checker::digest(field.to_string().as_bytes()),"field":field,"curve":curve,"sources":["docs/curves/standards/registry.json"]})
}
pub(super) fn catalog(inputs: &[Value]) -> R<(Vec<Value>, Vec<Value>)> {
    let mut rows: BTreeMap<String, Value> = BTreeMap::new();
    let mut findings = vec![];
    let mut seen = HashSet::new();
    for source in inputs {
        let key = source["source_id"].as_str().ok_or("missing source_id")?;
        if !seen.insert(key) {
            return Err(format!("duplicate source_id {key}"));
        }
        match import(source){
            Ok(row)=>{
                findings.push(json!({"source_id":key,"name":source["name"],"category":source["category"],"status":"imported","icv1_slug":row["slug"],
                    "ec1_status":if row["representations"].as_array().unwrap().is_empty(){"missing_full_generator"}else{"linked"},"source_record_sha256":checker::digest(source.to_string().as_bytes())}));
                let id=row["model_json"].as_str().unwrap().to_string();
                if let Some(existing)=rows.get_mut(&id){
                    if existing["order"]!=row["order"]{return Err(format!("conflicting orders for {key}"));}
                    for k in ["aliases","standard_names","standards_provenance","representations"] {for item in row[k].as_array().unwrap(){if !existing[k].as_array().unwrap().contains(item){existing[k].as_array_mut().unwrap().push(item.clone());}}}
                }else{rows.insert(id,row);}
            },
            Err(reason)=>findings.push(json!({"source_id":key,"name":source["name"],"category":source["category"],"status":"unresolved","reason":reason,"exists":null,"source_record_sha256":checker::digest(source.to_string().as_bytes())})),
        }
    }
    Ok((rows.into_values().collect(), findings))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn source() -> Value {
        json!({"source_id":"fixture/one","name":"fixture","category":"other","field":{"type":"Prime","p":"5"},"form":"Weierstrass","params":{"a":{"raw":"0"},"b":{"raw":"1"}},"order":"3","cofactor":"2","generator":{"x":{"raw":"0"},"y":{"raw":"1"}}})
    }
    #[test]
    fn imports_full_representation_and_accounts_for_gaps() {
        let first = source();
        let mut other = first.clone();
        other["source_id"] = json!("fixture/two");
        let mut gap = first.clone();
        gap["source_id"] = json!("fixture/gap");
        gap["field"]["type"] = json!("Extension");
        let (rows, findings) = catalog(&[first.clone(), other, gap]).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(findings.len(), 3);
        assert_eq!(findings[2]["status"], "unresolved");
        assert!(findings[2]["exists"].is_null());
        assert_eq!(rows[0]["representations"].as_array().unwrap().len(), 1);
        assert!(catalog(&[first.clone(), first]).is_err());
    }
    #[test]
    fn signed_coefficients_and_exact_integer_hashing() {
        let mut s = source();
        s["params"]["a"]["raw"] = json!("-1");
        assert_eq!(
            prime_parameter(&s, "a", &BigUint::from(5u32)).unwrap(),
            BigUint::from(4u32)
        );
        let integer = BigUint::one() << 300;
        assert_eq!(num(&integer).to_string(), integer.to_string());
        let mut bad = source();
        bad["generator"]["y"]["raw"] = json!("2");
        assert!(import(&bad).is_err());
    }
}
