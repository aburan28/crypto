//! Scoped catalogue and immutable synthetic-model index; no solving stage.
use crypto_lib::hash::sha256::sha256;
use num_bigint::{BigInt, BigUint};
use num_traits::One;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};
mod presentation;

const STUDY: &str = "research/iso1_weak_classes_20261007/two_branch_20261009";
const CATALOG: &str = "docs/curves/weak-families";
const OUTPUT: &str = "research/weak_curve_catalog_20261010";
const INPUTS: [(&str, &str); 5] = [
    ("research/iso1_weak_classes_20261007/two_branch_20261009/curve_records.json", "d16fb940ae8527ea5b44cad443970b1942bb943d9373f49d6ebab1479e15c790"),
    ("research/iso1_weak_classes_20261007/two_branch_20261009/p7_run1/stdout.txt", "fcf076d145a1b087c7dc4df26a6e3b32b39664bfe9350698956e68bfa607a786"),
    ("research/iso1_weak_classes_20261007/two_branch_20261009/p7_classes.csv", "ad6fcfe4422fa6961523f287ca3b56e5f1eb16976c157d89bb0dbfb3938dc595"),
    ("research/iso1_weak_classes_20261007/two_branch_20261009/summary.json", "e89831369d99ba74ab548e65a73a538a758d2f99f9aedd6ca36c381f7a21e789"),
    ("research/iso1_weak_classes_20261007/large_population_20261009/population_run1/curve_records.json", "7ce616701d8ae6d63e4e7dc289ad3914d27da59e2287f7abe7cc096fd978dbeb"),
];
fn digest(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}
fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn kv(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|s| s.split_once('='))
        .collect()
}
fn yn(b: bool) -> &'static str {
    if b {
        "yes"
    } else {
        "no"
    }
}

fn validate_catalog(v: &Value) -> Result<(), String> {
    if v["schema"] != "weak-curve-family-catalog/v1" || v["global_completeness"] != "open" {
        return Err("catalogue cannot promote scoped rules to global completeness".into());
    }
    let sources = v["sources"].as_array().ok_or("missing sources")?;
    let held: BTreeSet<_> = sources.iter().map(|s| s["id"].as_str().unwrap()).collect();
    if held.len() != sources.len() {
        return Err("duplicate source identifiers".into());
    }
    let mut ids = BTreeSet::new();
    for f in v["families"].as_array().ok_or("missing families")? {
        for k in [
            "id",
            "name",
            "category",
            "scope",
            "evidence",
            "domain",
            "recognition",
            "find_members",
            "exhaustive_scope",
            "subgroup_gate",
            "cost_gate",
            "local_status",
        ] {
            if f[k].as_str().is_none_or(|s| s.is_empty()) {
                return Err(format!("missing {k}"));
            }
        }
        if !ids.insert(f["id"].as_str().unwrap()) {
            return Err("duplicate family identifier".into());
        }
        if f["category"] == "structural_lead" && f["evidence"] != "structural_only" {
            return Err("structural signals require a separate witness before promotion".into());
        }
        if f["properties"]
            .as_array()
            .is_none_or(|a| a.is_empty() || a.iter().any(|p| p.as_str().is_none_or(str::is_empty)))
        {
            return Err("family properties must be nonempty statements".into());
        }
        for s in f["sources"].as_array().ok_or("family needs sources")? {
            if !held.contains(s.as_str().unwrap()) {
                return Err("dangling family source".into());
            }
        }
    }
    for r in v["relations"].as_array().ok_or("missing relations")? {
        if !ids.contains(r["from"].as_str().unwrap()) || !ids.contains(r["to"].as_str().unwrap()) {
            return Err("dangling relation".into());
        }
        if r["condition"].as_str().is_none_or(|s| s.is_empty()) {
            return Err("relation lacks hypotheses".into());
        }
    }
    Ok(())
}

fn verify_model(c: &Value) -> Result<(Value, BigUint), String> {
    let text = c["model_json"].as_str().ok_or("missing exact model")?;
    let hash = digest(text.as_bytes());
    if c["model_sha256"].as_str().is_some_and(|h| h != hash) {
        return Err("model digest mismatch".into());
    }
    let model: Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
    let p: BigUint = model["p"]
        .as_str()
        .unwrap()
        .parse()
        .map_err(|_| "invalid characteristic")?;
    let degree: u32 = model["k"]
        .as_str()
        .unwrap()
        .parse()
        .map_err(|_| "invalid degree")?;
    if degree == 0 || p <= BigUint::from(3u32) {
        return Err("outside frozen source domain".into());
    }
    let q = p.pow(degree);
    let t: BigInt = c["trace"]
        .as_str()
        .unwrap()
        .parse()
        .map_err(|_| "invalid trace")?;
    let n: BigUint = c
        .get("order")
        .or_else(|| c.get("group_order"))
        .unwrap()
        .as_str()
        .unwrap()
        .parse()
        .map_err(|_| "invalid order")?;
    if BigInt::from(n.clone()) != BigInt::from(q.clone()) + BigInt::one() - &t {
        return Err("trace/cardinality mismatch".into());
    }
    if &t * &t > BigInt::from(q.clone() * 4u32) {
        return Err("Hasse bound mismatch".into());
    }
    let id = c["icv1"].as_str().unwrap();
    let parts: Vec<_> = id.split(':').collect();
    let j = c["j_coefficients"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect::<Vec<_>>()
        .join(",");
    if parts.len() != 9
        || parts[1] != model["field"].as_str().unwrap()
        || parts[2] != c["trace"].as_str().unwrap()
        || parts[3] != n.to_string()
        || parts[4] != j
        || parts[8] != &hash[..12]
    {
        return Err("canonical identity fields mismatch".into());
    }
    if model["a"].as_array().unwrap().len() != degree as usize
        || model["b"].as_array().unwrap().len() != degree as usize
        || model["modulus"].as_array().unwrap().len() != degree as usize
    {
        return Err("field coefficient length mismatch".into());
    }
    Ok((model, q))
}

fn index(repo: &Path, catalog: &Value) -> Value {
    validate_catalog(catalog).unwrap();
    for source in catalog["sources"].as_array().unwrap() {
        if let Some(path) = source["path"].as_str() {
            assert!(repo.join(path).is_file(), "missing local source {path}");
        }
    }
    let mut bindings = Vec::new();
    for (path, expected) in INPUTS {
        let bytes = fs::read(repo.join(path)).unwrap();
        assert_eq!(digest(&bytes), expected, "frozen input changed: {path}");
        bindings.push(json!({"path":path,"sha256":expected,"bytes":bytes.len()}));
    }
    let mut oracle = BTreeMap::new();
    for line in fs::read_to_string(repo.join(STUDY).join("p7_classes.csv"))
        .unwrap()
        .lines()
        .skip(1)
    {
        let c: Vec<_> = line.split(',').collect();
        let t: i64 = c[1].parse().unwrap();
        let a = c[2] == "1";
        let b = c[3] == "1";
        assert_eq!(c[4] == "1", a || b);
        oracle.insert(t, (a, b));
    }
    assert_eq!(oracle.len(), 294);
    let classes: Vec<_> = oracle
        .iter()
        .map(|(t, (a, b))| {
            json!({
                "field_cardinality":"117649", "characteristic":"7", "degree":6,
                "trace":t.to_string(), "order":(117650-t).to_string(),
                "domain":"ordinary trace t = 2 modulo 4",
                "cubic_support":a, "quadratic_support":b, "combined_support":*a||*b,
                "evidence":INPUTS[2].0
            })
        })
        .collect();
    let class_counts = json!({"eligible":classes.len(),
        "cubic":oracle.values().filter(|(a,_)|*a).count(),
        "quadratic":oracle.values().filter(|(_,b)|*b).count(),
        "both":oracle.values().filter(|(a,b)|*a&&*b).count(),
        "combined":oracle.values().filter(|(a,b)|*a||*b).count(),
        "combined_zero":oracle.values().filter(|(a,b)|!*a&&!*b).count()});
    assert_eq!(class_counts["combined"], 248);
    assert_eq!(class_counts["combined_zero"], 46);
    let mut direct: BTreeMap<String, u8> = BTreeMap::new();
    let census = fs::read_to_string(repo.join(STUDY).join("p7_run1/stdout.txt")).unwrap();
    let field = kv(census.lines().next().unwrap());
    let modulus: Value = serde_json::from_str(field["modulus"]).unwrap();
    for line in census
        .lines()
        .filter(|l| l.starts_with("PLUS|") || l.starts_with("MINUS|"))
    {
        let m = kv(line);
        let t: i64 = m["trace"].parse().unwrap();
        if t % 7 == 0 {
            continue;
        }
        let j: Value = serde_json::from_str(m["j"]).unwrap();
        let flag = if line.starts_with("PLUS|") { 1 } else { 2 };
        *direct
            .entry(serde_json::to_string(&j).unwrap())
            .or_default() |= flag;
    }
    assert!(
        direct.values().all(|v| *v != 3),
        "ordinary cubic/quadratic models have different 2-torsion patterns"
    );
    let source = read_json(&repo.join(INPUTS[0].0));
    let population = read_json(&repo.join(INPUTS[4].0));
    let raw_summary = read_json(&repo.join(INPUTS[3].0));
    let searches: BTreeMap<_, _> = raw_summary["large_rows"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| (r["source"]["icv1"].as_str().unwrap().to_string(), r))
        .collect();
    assert_eq!(searches.len(), 512);
    let mut curves: BTreeMap<String, Value> = BTreeMap::new();
    for c in source["records"].as_array().unwrap() {
        let (model, q) = verify_model(c).unwrap();
        let p = model["p"].as_str().unwrap();
        let degree: u32 = model["k"].as_str().unwrap().parse().unwrap();
        let mut labels = json!({"model_g3_cubic":"outside_domain","model_g3_quadratic":"outside_domain","class_g3_cubic":"outside_domain","class_g3_quadratic":"outside_domain","combined_class":"outside_domain","other_family_membership":"unclassified","subgroup_advantage":"unmeasured"});
        if p == "7" && degree == 6 {
            assert_eq!(model["modulus"], modulus);
            let t: i64 = c["trace"].as_str().unwrap().parse().unwrap();
            let (a, b) = oracle[&t];
            let f = direct
                .get(&serde_json::to_string(&c["j_coefficients"]).unwrap())
                .copied()
                .unwrap_or(0);
            assert!(f & 1 == 0 || a);
            assert!(f & 2 == 0 || b);
            labels["model_g3_cubic"] = json!(yn(f & 1 != 0));
            labels["model_g3_quadratic"] = json!(yn(f & 2 != 0));
            labels["class_g3_cubic"] = json!(yn(a));
            labels["class_g3_quadratic"] = json!(yn(b));
            labels["combined_class"] = json!(if a || b {
                "geometric_support_exact"
            } else {
                "exact_family_zero"
            });
        } else {
            labels["higher_family_class"] = json!("geometric_control_witness");
            labels["direct_higher_family"] = json!("unclassified_in_this_index");
        }
        let row = json!({"slug":c["slug"],"icv1":c["icv1"],"model_sha256":digest(c["model_json"].as_str().unwrap().as_bytes()),"model_json":c["model_json"],"trace":c["trace"],"order":c["order"],"j_coefficients":c["j_coefficients"],"field":{"characteristic":p,"degree":degree,"modulus_low_coefficients":model["modulus"],"cardinality":q.to_string(),"bit_length":q.bits()},"subgroup_order":null,"generator":null,"endomorphism_order_conductor":null,"labels":labels,"inventory":"postselected_replayed_route_or_control","evidence":{"record_file":INPUTS[0].0,"coordinate_replay":format!("{STUDY}/{}",c["raw_source"].as_str().unwrap()),"class_oracle":if p=="7"&&degree==6{json!(INPUTS[2].0)}else{Value::Null},"direct_model_oracle":if p=="7"&&degree==6{json!(INPUTS[1].0)}else{Value::Null}}});
        assert!(curves
            .insert(c["icv1"].as_str().unwrap().to_string(), row)
            .is_none());
    }
    for r in population["records"].as_array().unwrap() {
        let c = r["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["role"] == "source")
            .unwrap();
        let (model, q) = verify_model(c).unwrap();
        assert_eq!(q.to_string(), r["field_cardinality"].as_str().unwrap());
        let degree: u32 = model["k"].as_str().unwrap().parse().unwrap();
        let search = searches[c["icv1"].as_str().unwrap()];
        assert_eq!(search["source"]["trace"], c["trace"]);
        assert_eq!(search["source"]["order"], c["group_order"]);
        assert_ne!(search["cost"]["status"], "VERIFIED_WITNESS");
        let labels = json!({"model_g3_cubic":if degree==6{"unclassified_in_this_index"}else{"outside_domain"},"model_g3_quadratic":if degree==6{"unclassified_in_this_index"}else{"outside_domain"},"class_g3_cubic":if degree==6&&!r["admitted"].as_bool().unwrap(){"excluded_by_proved_branch_condition"}else if degree==6{"unresolved"}else{"outside_domain"},"class_g3_quadratic":if degree==6{"unresolved"}else{"outside_domain"},"combined_class":if degree==6{"unresolved"}else{"outside_domain"},"higher_family_class":if degree!=6{"unresolved"}else{"outside_domain"},"other_family_membership":"unclassified","subgroup_advantage":"unmeasured"});
        let row = json!({"slug":c["slug"],"icv1":c["icv1"],"model_sha256":digest(c["model_json"].as_str().unwrap().as_bytes()),"model_json":c["model_json"],"trace":c["trace"],"order":c["group_order"],"j_coefficients":c["j_coefficients"],"field":{"characteristic":model["p"],"degree":degree,"modulus_low_coefficients":model["modulus"],"cardinality":q.to_string(),"bit_length":q.bits()},"subgroup_order":null,"generator":null,"endomorphism_order_conductor":null,"labels":labels,"inventory":"independent_model_law_large_source","seed":r["seed"],"historical_cubic_admission":r["admitted"],"historical_cubic_label":r["class_label"],"bounded_search":search["cost"],"evidence":{"record_file":INPUTS[4].0,"bounded_search_summary":INPUTS[3].0,"raw_output":format!("research/iso1_weak_classes_20261007/large_population_20261009/population_run1/{}",r["raw_output"].as_str().unwrap())}});
        assert!(curves
            .insert(c["icv1"].as_str().unwrap().to_string(), row)
            .is_none());
    }
    assert_eq!(curves.len(), 814);
    let models: Vec<_> = curves.into_values().collect();
    let mut counts = BTreeMap::<String, usize>::new();
    for c in &models {
        *counts
            .entry(format!("inventory:{}", c["inventory"].as_str().unwrap()))
            .or_default() += 1;
        for key in [
            "model_g3_cubic",
            "model_g3_quadratic",
            "combined_class",
            "higher_family_class",
        ] {
            if let Some(status) = c["labels"][key].as_str() {
                *counts.entry(format!("{key}:{status}")).or_default() += 1;
            }
        }
    }
    let mut categories = BTreeMap::<String, usize>::new();
    for f in catalog["families"].as_array().unwrap() {
        *categories
            .entry(f["category"].as_str().unwrap().into())
            .or_default() += 1;
    }
    json!({"schema":"weak-curve-model-index/v1","family_catalog_sha256":digest(&fs::read(repo.join(CATALOG).join("catalog.json")).unwrap()),"input_bindings":bindings,"summary":{"families":catalog["families"].as_array().unwrap().len(),"models":models.len(),"categories":categories,"inventory_counts":counts,"class_oracle_counts":class_counts,"count_scope":"inventory counts; postselected replay models and independent sources are separate input laws","global_completeness":"open"},"class_oracle":classes,"models":models})
}

fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "build|check|render REPO");
    let repo = Path::new(&args[1]);
    let catalog = read_json(&repo.join(CATALOG).join("catalog.json"));
    let data = index(repo, &catalog);
    let text = serde_json::to_vec_pretty(&data).unwrap();
    match args[0].as_str() {
        "build" => {
            fs::write(repo.join(CATALOG).join("models.json"), text).unwrap();
            fs::write(
                repo.join(OUTPUT).join("summary.json"),
                serde_json::to_vec_pretty(&data["summary"]).unwrap(),
            )
            .unwrap();
        }
        "check" => {
            assert_eq!(
                fs::read(repo.join(CATALOG).join("models.json")).unwrap(),
                text
            );
            assert_eq!(
                read_json(&repo.join(OUTPUT).join("summary.json")),
                data["summary"]
            );
            presentation::check(repo, &catalog, &data);
        }
        "render" => presentation::render(repo, &catalog, &data),
        _ => panic!("unknown catalogue operation"),
    }
    println!(
        "CATALOGUE_VERIFIED|families={}|models={}|large_class_labels=unresolved",
        data["summary"]["families"], data["summary"]["models"]
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    fn catalogue() -> Value {
        read_json(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../../docs/curves/weak-families/catalog.json"),
        )
    }
    #[test]
    fn global_completeness_promotion_is_rejected() {
        let mut c = catalogue();
        assert!(validate_catalog(&c).is_ok());
        c["global_completeness"] = json!("complete");
        assert!(validate_catalog(&c).is_err());
    }
    #[test]
    fn structural_property_needs_a_separate_witness() {
        let mut c = catalogue();
        let f = c["families"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|f| f["category"] == "structural_lead")
            .unwrap();
        f["evidence"] = json!("native_verified_geometry");
        assert!(validate_catalog(&c).is_err());
    }
    #[test]
    fn trace_and_model_mutations_fail_closed() {
        let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
        let mut c = read_json(&repo.join(INPUTS[0].0))["records"][0].clone();
        assert!(verify_model(&c).is_ok());
        c["order"] = json!("117381");
        assert!(verify_model(&c).is_err());
        let mut c = read_json(&repo.join(INPUTS[0].0))["records"][0].clone();
        c["model_json"] = json!(c["model_json"]
            .as_str()
            .unwrap()
            .replace("\"v\":\"1\"", "\"v\":\"2\""));
        assert!(verify_model(&c).is_err());
    }
}
