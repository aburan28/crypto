//! Complete pinned catalogue with explicit operation capabilities.
use isogeny_algos::{api::PrimeCurve, bigint::Big, int::Int};
use serde_json::{json, Value};
use std::sync::OnceLock;

pub fn registry() -> &'static Value {
    static DATA: OnceLock<Value> = OnceLock::new();
    DATA.get_or_init(|| serde_json::from_str(crate::REGISTRY).unwrap())
}
pub fn integer(value: &Value) -> Int {
    let value = value
        .as_str()
        .map(str::to_owned)
        .unwrap_or_else(|| value.to_string());
    let (negative, digits) = value
        .strip_prefix('-')
        .map_or((false, value.as_str()), |v| (true, v));
    let number = if let Some(h) = digits.strip_prefix("0x") {
        Big::from_dec(
            &num_bigint::BigUint::parse_bytes(h.as_bytes(), 16)
                .unwrap()
                .to_str_radix(10),
        )
    } else {
        Big::from_dec(digits)
    };
    let number = Int::from_big(&number);
    if negative {
        -&number
    } else {
        number
    }
}
fn normalized(name: &str) -> String {
    name.chars()
        .filter(char::is_ascii_alphanumeric)
        .map(|c| c.to_ascii_lowercase())
        .collect()
}
pub fn resolve(name: &str) -> Result<&'static Value, String> {
    let rows = registry()["curves"].as_array().unwrap();
    if let Some(row) = rows.iter().find(|r| r["slug"] == name) {
        return Ok(row);
    }
    let needle = normalized(name);
    let matches: Vec<_> = rows
        .iter()
        .filter(|r| {
            ["aliases", "standard_names"].iter().any(|key| {
                r[*key].as_array().is_some_and(|names| {
                    names
                        .iter()
                        .any(|n| n.as_str().is_some_and(|s| normalized(s) == needle))
                })
            })
        })
        .collect();
    match matches.as_slice() {
        [row] => Ok(*row),
        [] => Err(format!("unknown curve {name}; run isogeny curves")),
        _ => Err(format!(
            "ambiguous curve {name}; select its exact ICV1 slug"
        )),
    }
}
pub fn label(row: &'static Value) -> &'static str {
    let names = row["standard_names"].as_array().unwrap();
    if names.iter().any(|n| n == "P-192") {
        "p192"
    } else if names.iter().any(|n| n == "P-224") {
        "p224"
    } else {
        row["slug"].as_str().unwrap()
    }
}
pub fn field_order(row: &Value) -> Result<Int, String> {
    let p = &row["params"];
    match row["family"].as_str().unwrap_or("") {
        "prime" => Ok(integer(&p["p"])),
        "binary" | "subfield" | "koblitz" => {
            let bits = p["m"]
                .as_u64()
                .or_else(|| p["n"].as_u64())
                .ok_or("missing binary field degree")?;
            Ok(Int::from(2i64).pow(bits.try_into().map_err(|_| "field degree is too large")?))
        }
        "extension" => {
            let base = integer(&p["p"]);
            let k = p["k"].as_u64().ok_or("missing extension degree")?;
            Ok(base.pow(k.try_into().map_err(|_| "field degree is too large")?))
        }
        _ => Err("field adapter is not available".into()),
    }
}
pub fn source_order(name: &str) -> Int {
    integer(&resolve(name).unwrap()["order"])
}
pub fn prime_curve(row: &Value) -> Result<PrimeCurve, String> {
    if row["family"] != "prime" {
        return Err("construction adapter is not yet implemented for this field family".into());
    }
    let p = &row["params"];
    if p["a"].is_null() || p["b"].is_null() {
        return Err(
            "this model requires a checked coordinate conversion before construction".into(),
        );
    }
    PrimeCurve::new(integer(&p["p"]), integer(&p["a"]), integer(&p["b"]))
}
pub fn representation(row: &Value) -> Option<&Value> {
    row["representations"].as_array()?.iter().find(|r| {
        r["curve"]["model"] == "short Weierstrass"
            && r["curve"]["generator"]
                .as_array()
                .is_some_and(|g| g.len() == 2 && g.iter().all(|v| !v.is_null()))
            && !r["curve"]["subgroup_order"].is_null()
    })
}
pub fn capabilities(row: &Value) -> Value {
    let bits = field_order(row).map(|q| q.bits()).ok();
    let prime = row["family"] == "prime";
    let short = prime && !row["params"]["a"].is_null() && !row["params"]["b"].is_null();
    let construction = short && bits.is_some_and(|n| n <= 640);
    let replay = construction && representation(row).is_some();
    let reason = if !prime {
        "construction requires a binary or extension-field map adapter"
    } else if !short {
        "construction requires checked Montgomery/Edwards coordinate conversion"
    } else if !construction {
        "field width exceeds the constructor limit"
    } else if !replay {
        "no complete registered public subgroup generator"
    } else {
        "registered short-Weierstrass model and public subgroup representation"
    };
    json!({"screen":bits.is_some(),"construction":construction,"independent_replay":replay,"field_bits":bits,"detail":reason})
}
pub fn inventory() -> Value {
    let rows: Vec<_> = registry()["curves"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            json!({
        "icv1":row["icv1"],"slug":row["slug"],"family":row["family"],
        "aliases":row["aliases"],"standard_names":row["standard_names"],
        "model_form":row["params"]["model_form"],"sources":row["sources"],
        "capabilities":capabilities(row),"representations":row["representations"]})
        })
        .collect();
    json!({"schema":"isogeny-catalogue/v1","curves":rows,"count":rows.len(),
        "standards_import_coverage":serde_json::from_str::<Value>(include_str!("../../../docs/curves/standards/coverage.json")).unwrap(),
        "scope":"all models in the pinned repository inventory; not a worldwide completeness claim"})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_models_are_accounted_for_with_operation_capabilities() {
        let data = inventory();
        assert_eq!(data["count"], 345);
        for row in data["curves"].as_array().unwrap() {
            assert_eq!(
                resolve(row["slug"].as_str().unwrap()).unwrap()["slug"],
                row["slug"]
            );
            assert!(row["capabilities"]["screen"].is_boolean());
            assert!(row["capabilities"]["construction"].is_boolean());
        }
    }
    #[test]
    fn common_standard_names_resolve_without_merging_different_models() {
        for name in [
            "P-192",
            "P-224",
            "P-256",
            "P-384",
            "P-521",
            "secp256k1",
            "brainpoolP384r1",
            "SM2",
            "Curve25519",
            "Curve448",
            "K-163",
        ] {
            assert!(resolve(name).is_ok(), "{name}");
        }
        assert_eq!(label(resolve("P-192").unwrap()), "p192");
        assert!(capabilities(resolve("P-521").unwrap())["construction"]
            .as_bool()
            .unwrap());
        assert_eq!(
            capabilities(resolve("K-163").unwrap())["construction"],
            false
        );
    }
}
