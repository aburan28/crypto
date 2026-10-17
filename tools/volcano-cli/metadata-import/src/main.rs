//! Metadata-only import of native plan identities. No curve arithmetic or counting.
//! Carries existing arbitrary-precision JSON integers and property order unchanged.
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

fn write(path: &Path, value: &Value, indent: &[u8]) {
    let mut out = vec![];
    let mut serializer = serde_json::Serializer::with_formatter(
        &mut out,
        serde_json::ser::PrettyFormatter::with_indent(indent),
    );
    serde::Serialize::serialize(value, &mut serializer).unwrap();
    out.push(b'\n');
    fs::write(path, out).unwrap();
}
fn normalized(s: &str) -> String {
    s.trim()
        .trim_matches('`')
        .chars()
        .filter(|c| !c.is_whitespace() && !matches!(c, '{' | '}' | '_'))
        .flat_map(char::to_lowercase)
        .collect()
}
fn collect_json(path: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(path).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            collect_json(&path, out);
        } else if path.extension().is_some_and(|s| s == "json") {
            out.push(path);
        }
    }
}
fn records<'a>(v: &'a Value, out: &mut Vec<&'a Value>) {
    if v.get("p").is_some()
        && v.get("a").is_some()
        && v.get("b").is_some()
        && v.get("order").is_some_and(|x| !x.is_null())
    {
        out.push(v);
    }
    match v {
        Value::Array(xs) => {
            for x in xs {
                records(x, out)
            }
        }
        Value::Object(xs) => {
            for (_, x) in xs {
                records(x, out)
            }
        }
        _ => {}
    }
}
fn integer(v: &Value) -> String {
    v.as_str()
        .map(str::to_owned)
        .unwrap_or_else(|| v.to_string())
}
fn number_key(s: &str) -> (usize, &str) {
    (s.len(), s)
}
fn rank(v: &Value) -> (u8, u64, String, u64, String) {
    let family = v["family"].as_str().unwrap();
    let params = &v["params"];
    let r = match family {
        "koblitz" => 0,
        "subfield" => 1,
        "binary" => 2,
        "prime" => 3,
        "extension" => 4,
        _ => panic!("unknown registry family"),
    };
    (
        r,
        params
            .get("n")
            .or_else(|| params.get("m"))
            .and_then(Value::as_u64)
            .unwrap_or(0),
        params["p"].as_str().unwrap_or("0").into(),
        if r == 0 {
            params["a"].as_u64().unwrap_or(0)
        } else {
            0
        },
        v["slug"].as_str().unwrap().into(),
    )
}
fn main() {
    let repo = PathBuf::from(std::env::args().nth(1).expect("repository root"));
    let study = repo.join(
        std::env::args()
            .nth(2)
            .unwrap_or_else(|| "research/full_volcano_classes_20261016".into()),
    );
    let registry_path = repo.join("docs/curves/registry.json");
    let mut registry: Value = serde_json::from_slice(&fs::read(&registry_path).unwrap()).unwrap();
    let mut files = vec![];
    collect_json(&study, &mut files);
    files.sort();
    let docs = files
        .iter()
        .map(|p| {
            (
                p,
                serde_json::from_slice::<Value>(&fs::read(p).unwrap()).unwrap(),
            )
        })
        .collect::<Vec<_>>();
    let mut added = vec![];
    for (path, plan) in &docs {
        if plan["schema"] != "ordinary-isogeny-class/v1"
            || plan["source_identity"]["model_json"].is_null()
        {
            continue;
        }
        let identity = &plan["source_identity"];
        let model: Value = serde_json::from_str(identity["model_json"].as_str().unwrap()).unwrap();
        let id = identity["icv1"].as_str().unwrap();
        let fields = id.split(':').collect::<Vec<_>>();
        assert_eq!(fields.len(), 9);
        let p = integer(&model["p"]);
        let a = integer(&model["a"]);
        let b = integer(&model["b"]);
        let order = fields[3];
        let mut sources = vec![];
        for (file, doc) in &docs {
            let mut rs = vec![];
            records(doc, &mut rs);
            if rs.iter().any(|r| {
                integer(&r["p"]) == p
                    && integer(&r["a"]) == a
                    && integer(&r["b"]) == b
                    && integer(&r["order"]) == order
            }) {
                sources.push(
                    file.strip_prefix(&repo)
                        .unwrap()
                        .to_string_lossy()
                        .to_string(),
                );
            }
        }
        if sources.is_empty() {
            sources.push(
                path.strip_prefix(&repo)
                    .unwrap()
                    .to_string_lossy()
                    .to_string(),
            );
        }
        sources.truncate(6);
        let curves = registry["curves"].as_array_mut().unwrap();
        if let Some(existing) = curves
            .iter_mut()
            .find(|c| c["model_json"] == identity["model_json"])
        {
            assert_eq!(existing["icv1"], identity["icv1"]);
            existing["sources"] = serde_json::to_value(&sources).unwrap();
            continue;
        }
        curves.push(json!({"slug":identity["slug"],"icv1":identity["icv1"],"family":"prime",
            "params":{"p":p,"a":a,"b":b},"model_json":identity["model_json"],"trace":fields[2].parse::<i64>().unwrap(),
            "order":order,"j":fields[4],"end":"unk","aliases":[],"standard_names":[],"sources":sources,
            "representations":[],"ec1_unresolved":"no tracked record carries this model's subgroup and generator"}));
        added.push(identity["slug"].clone());
    }
    registry["curves"].as_array_mut().unwrap().sort_by(|a, b| {
        let a = rank(a);
        let b = rank(b);
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
            .then(number_key(&a.2).cmp(&number_key(&b.2)))
            .then(a.3.cmp(&b.3))
            .then(a.4.cmp(&b.4))
    });
    write(&registry_path, &registry, b" ");
    let alias_path = repo.join("src/cryptanalysis/curve_aliases.json");
    let mut aliases: Value = serde_json::from_slice(&fs::read(&alias_path).unwrap()).unwrap();
    let mut map: BTreeMap<String, Value> =
        serde_json::from_value(aliases["aliases"].clone()).unwrap();
    for c in registry["curves"].as_array().unwrap() {
        for key in ["icv1", "slug"] {
            map.insert(normalized(c[key].as_str().unwrap()), c["slug"].clone());
        }
    }
    aliases["aliases"] = serde_json::to_value(map).unwrap();
    write(&alias_path, &aliases, b"");
    println!(
        "{}",
        json!({"status":"METADATA_IMPORTED","added":added,"source":"frozen native plans; no arithmetic or subgroup inferred"})
    );
}
