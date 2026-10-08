//! Content-addressed geometric cover graph. EC1 bindings identify applicable
//! representations; they make no assertion about Jacobian subgroup transport.
use super::checker::{digest, number};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};
use std::collections::BTreeMap;

fn integer(v: &Value) -> Result<BigUint, String> {
    number(
        &v.as_str()
            .map(str::to_owned)
            .unwrap_or_else(|| v.to_string()),
    )
}

fn bind(model: &Value, rep: &Value, row: &Value) -> Result<(), String> {
    let field = &rep["field"];
    let curve = &rep["curve"];
    let c = &curve["coefficients"];
    let eq = |a: &Value, b: &Value| -> Result<bool, String> { Ok(integer(a)? == integer(b)?) };
    let valid = match model["form"].as_str() {
        Some("y^2+xy=x^3+a*x^2+b") => {
            let mut modulus = BigUint::zero();
            for exp in field["modulus_exponents"]
                .as_array()
                .ok_or("missing EC1 basis")?
            {
                let power = exp.as_u64().ok_or("invalid basis exponent")?;
                if power > 4096 || modulus.bit(power) {
                    return Err("invalid EC1 basis".into());
                }
                modulus.set_bit(power, true);
            }
            field["characteristic"] == 2
                && field["representation"] == "polynomial"
                && field["degree"].as_u64() == Some(modulus.bits() - 1)
                && modulus == integer(&model["modulus"])?
                && curve["model"] == "binary Weierstrass"
                && c.as_array().is_some_and(|a| a.len() == 5)
                && integer(&c[0])?.is_one()
                && eq(&c[1], &model["a"])?
                && integer(&c[2])?.is_zero()
                && integer(&c[3])?.is_zero()
                && eq(&c[4], &model["b"])?
        }
        Some(form) => {
            if field["representation"] != "prime"
                || field["degree"] != 1
                || !eq(&field["characteristic"], &model["p"])?
            {
                return Err("EC1 field differs from model".into());
            }
            match form {
                "y^2=x^3+a*x+b" => {
                    curve["model"] == "short Weierstrass"
                        && c.as_array().is_some_and(|a| a.len() == 5)
                        && integer(&c[0])?.is_zero()
                        && integer(&c[1])?.is_zero()
                        && integer(&c[2])?.is_zero()
                        && eq(&c[3], &model["a"])?
                        && eq(&c[4], &model["b"])?
                }
                "B*y^2=x^3+A*x^2+x" => {
                    curve["model"] == "Montgomery"
                        && eq(&c["a"], &model["A"])?
                        && eq(&c["b"], &model["B"])?
                }
                "a*x^2+y^2=1+d*x^2*y^2" => {
                    curve["model"] == "TwistedEdwards"
                        && eq(&c["a"], &model["a"])?
                        && eq(&c["d"], &model["d"])?
                }
                "x^2+y^2=c^2*(1+d*x^2*y^2)" => {
                    curve["model"] == "Edwards"
                        && eq(&c["c"], &model["c"])?
                        && eq(&c["d"], &model["d"])?
                }
                _ => false,
            }
        }
        _ => false,
    };
    if !valid
        || integer(&curve["subgroup_order"])? * integer(&curve["cofactor"])?
            != integer(&row["order"])?
    {
        return Err("EC1 tuple does not describe the linked ICV1 model/order".into());
    }
    Ok(())
}

pub(super) fn graph(registry: &Value, report: &Value) -> Result<Value, String> {
    let mut covers = BTreeMap::new();
    let mut maps = BTreeMap::new();
    let mut curves = vec![];
    let findings = report["curves"].as_array().ok_or("missing findings")?;
    let rows = registry["curves"]
        .as_array()
        .ok_or("missing registry curves")?;
    if rows.len() != findings.len() {
        return Err("registry/finding cardinality mismatch".into());
    }
    for (row, finding) in rows.iter().zip(findings) {
        let raw = row["model_json"].as_str().ok_or("missing model")?;
        let model: Value = serde_json::from_str(raw).map_err(|e| e.to_string())?;
        let model_hash = digest(raw.as_bytes());
        if finding["slug"] != row["slug"] || finding["model_sha256"] != model_hash {
            return Err("cover join mismatch".into());
        }
        let mut representations = vec![];
        if let Some(reps) = row["representations"].as_array() {
            for rep in reps {
                bind(&model, rep, row)?;
                let preimage = json!({"field":rep["field"],"curve":rep["curve"]});
                let hash = digest(preimage.to_string().as_bytes());
                let uid = format!("urn:ec-record:1:sha256:{hash}");
                if rep["curve_uid"] != uid
                    || !rep["ec1"]
                        .as_str()
                        .is_some_and(|s| s.ends_with(&format!("h{}", &hash[..12])))
                {
                    return Err(format!("EC1 preimage mismatch for {}", row["slug"]));
                }
                representations.push(json!({"ec1":rep["ec1"],"curve_uid":uid,
                    "field":rep["field"],"curve":rep["curve"],"subgroup_transport":{"status":"not_tested","evidence_ref":null}}));
            }
        }
        let mut refs = vec![];
        if finding["exists"] == true {
            let c = &finding["certificate"];
            let mut field = json!({"label":model["field"]});
            for k in ["p", "modulus"] {
                if !model[k].is_null() {
                    field[k] = model[k].clone();
                }
            }
            let preimage = json!({"schema_version":"hyperelliptic-model/v1","field":field,"form":"v^2+h(u)*v=f(u)","h":c["h"],"f":c["f"],"genus":c["genus"]});
            let hash = digest(preimage.to_string().as_bytes());
            let hc_uid = format!("urn:hc-model:1:sha256:{hash}");
            covers.insert(hc_uid.clone(),json!({"cover_uid":hc_uid,"cover_id":format!("HC1G{}h{}",c["genus"],&hash[..12]),"record":preimage}));
            let map_record = json!({"schema_version":"curve-cover-map/v1","source_cover_uid":hc_uid,
                "target_model_sha256":model_hash,"target_model_json":raw,"degree":c["degree"],
                "x":c["x"],"y_v":c["y_v"],"y_0":c["y_0"],"target_model_map":c["target_model_map"]});
            let map_hash = digest(map_record.to_string().as_bytes());
            let map_uid = format!("urn:curve-cover-map:1:sha256:{map_hash}");
            maps.insert(map_uid.clone(),json!({"map_uid":map_uid,"map_id":format!("CV1D{}h{}",c["degree"],&map_hash[..12]),
                "record":map_record,"verification_status":finding["status"],"field_check":finding["field_check"],"proof_refs":["docs/curves/COVERS.md"]}));
            refs.push(json!({"cover_uid":hc_uid,"map_uid":map_uid,"direction":"H -> E"}));
        }
        curves.push(json!({"icv1_identity":{"slug":row["slug"],"full":row["icv1"]},"model_sha256":model_hash,
            "model_json":raw,"representations":representations,"cover_links":refs,"status":finding["status"],
            "standard_names":row["standard_names"],"standards_provenance":row["standards_provenance"]}));
    }
    Ok(
        json!({"schema_version":"curve-cover-graph/v1","identity_rule":"SHA-256 of sorted-key compact UTF-8 JSON of each mathematical record; metadata excluded",
        "registry_sha256":report["registry_sha256"],"global_standard_coverage":"not_claimed; see standards/coverage.json",
        "curves":curves,"covers":covers,"maps":maps}),
    )
}
