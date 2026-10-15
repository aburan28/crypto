//! Incremental native index writer. Existing raw rows remain byte-for-byte
//! intact unless this study adds a source/representation to that model.
use super::*;
use crypto_lib::cryptanalysis::{curve_id, ecbench::canonical::digest};
use num_bigint::BigUint;
fn q(v: &str) -> String {
    serde_json::to_string(v).unwrap()
}
fn object_keys(v: &Value) -> Vec<String> {
    let m = v.as_object().unwrap();
    let order: &[&str] = if m.contains_key("slug") {
        &[
            "slug",
            "icv1",
            "family",
            "params",
            "model_json",
            "trace",
            "order",
            "j",
            "end",
            "aliases",
            "standard_names",
            "sources",
            "representations",
            "ec1_unresolved",
        ]
    } else if m.contains_key("curve_uid") {
        &[
            "ec1",
            "curve_uid",
            "field_sha256",
            "field",
            "curve",
            "sources",
        ]
    } else if m.contains_key("representation") {
        &[
            "characteristic",
            "degree",
            "representation",
            "element_encoding",
        ]
    } else if m.contains_key("coefficients") {
        &[
            "model",
            "coefficients",
            "subgroup_order",
            "cofactor",
            "generator",
            "target_group",
        ]
    } else if m.contains_key("p") && m.contains_key("a") {
        &["p", "a", "b", "generator_call"]
    } else {
        &[]
    };
    let mut keys: Vec<String> = order
        .iter()
        .filter(|k| m.contains_key(**k))
        .map(|x| x.to_string())
        .collect();
    for k in m.keys() {
        if !keys.contains(k) {
            keys.push(k.clone());
        }
    }
    keys
}
fn render(v: &Value, depth: usize) -> String {
    let pad = " ".repeat(depth);
    let inner = " ".repeat(depth + 1);
    match v {
        Value::Object(m) => {
            let x: Vec<_> = object_keys(v)
                .iter()
                .map(|k| format!("{inner}{}: {}", q(k), render(&m[k], depth + 1)))
                .collect();
            if x.is_empty() {
                "{}".into()
            } else {
                format!("{{\n{}\n{pad}}}", x.join(",\n"))
            }
        }
        Value::Array(a) => {
            if a.is_empty() {
                "[]".into()
            } else {
                format!(
                    "[\n{}\n{pad}]",
                    a.iter()
                        .map(|x| format!("{inner}{}", render(x, depth + 1)))
                        .collect::<Vec<_>>()
                        .join(",\n")
                )
            }
        }
        _ => v.to_string(),
    }
}
fn segments(text: &str, key: &str) -> (usize, usize, Vec<String>) {
    let pos = text.find(&format!("\"{key}\": [")).unwrap() + key.len() + 4;
    let bytes = text.as_bytes();
    assert_eq!(bytes[pos], b'[');
    let mut quote = false;
    let mut escape = false;
    let mut depth = 0;
    let mut start = 0;
    let mut rows = Vec::new();
    for i in pos + 1..bytes.len() {
        let b = bytes[i];
        if quote {
            if escape {
                escape = false
            } else if b == b'\\' {
                escape = true
            } else if b == b'"' {
                quote = false
            }
            continue;
        }
        match b {
            b'"' => quote = true,
            b'{' => {
                if depth == 0 {
                    start = i;
                }
                depth += 1;
            }
            b'}' => {
                depth -= 1;
                if depth == 0 {
                    rows.push(text[start..=i].to_string());
                }
            }
            b']' if depth == 0 => return (pos, i, rows),
            _ => (),
        }
    }
    panic!("unterminated index array")
}
fn big(v: &Value) -> BigUint {
    let t = v
        .as_str()
        .map(str::to_string)
        .unwrap_or_else(|| v.to_string());
    if let Some(h) = t.strip_prefix("0x") {
        BigUint::parse_bytes(h.as_bytes(), 16).unwrap()
    } else {
        t.parse().unwrap()
    }
}

pub fn dashboard() {
    let tp = "docs/ic/progress-timeline.json";
    let mut timeline = read(tp);
    timeline["updated"] = json!("2026-10-15");
    let mut pending = timeline["non_numeric_results"]
        .as_array()
        .cloned()
        .unwrap_or_default();
    pending.retain(|v| v["anchor"] != "#conductor-gap-native-20261015");
    pending.push(json!({"date":"2026-10-15","anchor":"#conductor-gap-native-20261015","label":"Conductor-gap native orbit and full-solve controls","value":null,"status":"All 500 primary and 265 replication attempts verified; L0; native work incompletely priced","source":"https://github.com/aburan28/crypto/blob/main/research/conductor_gap_cm_orbits_20261015/REPORT.md","reason_not_plotted":"No complete cold one-target IC/rho ratio exists. Zero charged GAE is an incomplete lower bound, not a measured zero cost."}));
    timeline["non_numeric_results"] = json!(pending);
    write(tp, &timeline);
    let path = "docs/index-calculus-scoreboard.html";
    let mut html = fs::read_to_string(path).unwrap();
    let panel = include_str!("dashboard_panel.html");
    if let Some(start) = html.find("<!-- conductor-gap-start -->") {
        let end = html[start..].find("<!-- conductor-gap-end -->").unwrap()
            + start
            + "<!-- conductor-gap-end -->".len();
        html.replace_range(start..end, panel.trim_end());
    } else {
        let before = "<section class=\"dash-section\" id=\"shared-rank-folded-table-20261003\"";
        let index = html.find(before).unwrap();
        html.insert_str(index, panel);
    }
    let mark = "<script type=\"application/json\" id=\"progress-data\">";
    let start = html.find(mark).unwrap() + mark.len();
    let end = html[start..].find("</script>").unwrap() + start;
    html.replace_range(start..end, fs::read_to_string(tp).unwrap().trim_end());
    if !html.contains("conductor-gap-unplotted") {
        let start = html.find(mark).unwrap();
        html.insert_str(start,"<p id=\"conductor-gap-unplotted\">15 October: <a href=\"#conductor-gap-native-20261015\">conductor-gap native controls</a> have no numerical IC/rho point while native costs remain incompletely priced. The canonical timeline retains this result in its non-numeric ledger.</p>\n");
    }
    fs::write(path, html).unwrap();
}
fn integer(s: &str) -> Value {
    Value::Number(s.parse().unwrap())
}
fn scan(v: &Value, source: &str, rows: &mut BTreeMap<String, (Value, Vec<String>)>) {
    if let Some(m) = v.as_object() {
        if m.contains_key("p")
            && m.contains_key("a")
            && m.contains_key("b")
            && (m.contains_key("group_order") || m.contains_key("order"))
        {
            let p = big(&v["p"]);
            let a = big(&v["a"]) % &p;
            let b = big(&v["b"]) % &p;
            let n = big(v.get("group_order").unwrap_or(&v["order"]));
            let id = curve_id::prime(&p, &a, &b, &n).unwrap();
            let trace: num_bigint::BigInt =
                num_bigint::BigInt::from(p.clone()) + 1 - num_bigint::BigInt::from(n.clone());
            let fields: Vec<_> = id.icv1.split(':').collect();
            let entry = json!({"slug":id.slug,"icv1":id.icv1,"family":"prime","params":{"p":p.to_string(),"a":a.to_string(),"b":b.to_string()},"model_json":id.model_json,"trace":integer(&trace.to_string()),"order":n.to_string(),"j":fields[4],"end":fields[5],"aliases":[],"standard_names":[],"sources":[],"representations":[]});
            let e = rows.entry(id.slug.clone()).or_insert((entry, Vec::new()));
            assert_eq!(e.0["icv1"], id.icv1, "short display slug collision");
            assert_eq!(e.0["model_json"], id.model_json, "model identity collision");
            if !e.1.contains(&source.to_string()) && e.1.len() < 6 {
                e.1.push(source.into());
            }
            if let Some(g) = v.get("generator").and_then(|g| g.as_array()) {
                let r = big(&v["r"]);
                let field = json!({"characteristic":integer(&p.to_string()),"degree":1,"representation":"prime","element_encoding":"hex integer modulo the characteristic"});
                let curve = json!({"model":"short Weierstrass","coefficients":[0,0,0,integer(&a.to_string()),integer(&b.to_string())],"subgroup_order":r.to_string(),"cofactor":integer(&(&n/&r).to_string()),"generator":g,"target_group":"prime-order subgroup"});
                let h = digest(&json!({"field":field,"curve":curve})).unwrap();
                let rep = json!({"ec1":format!("EC1P{}Cfph{}",p.bits(),&h[..12]),"curve_uid":format!("urn:ec-record:1:sha256:{h}"),"field_sha256":digest(&field).unwrap(),"field":field,"curve":curve,"sources":[source]});
                let reps = e.0["representations"].as_array_mut().unwrap();
                if !reps.iter().any(|x| x["curve_uid"] == rep["curve_uid"]) {
                    reps.push(rep);
                }
            }
        }
        for x in m.values() {
            scan(x, source, rows)
        }
    } else if let Some(a) = v.as_array() {
        for x in a {
            scan(x, source, rows)
        }
    }
}
pub fn run() {
    let root = Path::new(ROOT);
    let catalog = read(root.join("prepared/catalog.json"));
    let records:Vec<_>=catalog.as_array().unwrap().iter().map(|v|{let s=&v["spec"];json!({"p":s["p"],"a":s["a"],"b":s["b"],"order":s["group_order"],"r":s["r"],"generator":[format!("0x{:x}",n(s,"gx")),format!("0x{:x}",n(s,"gy"))]})}).collect();
    write(root.join("prepared/curve_ids.json"), &json!(records));
    let listing = std::process::Command::new("git")
        .args(["ls-files", "--", ROOT])
        .output()
        .unwrap();
    assert!(listing.status.success());
    let mut files: Vec<String> = String::from_utf8(listing.stdout)
        .unwrap()
        .lines()
        .filter(|p| p.ends_with(".json"))
        .map(str::to_string)
        .collect();
    files.sort();
    let mut additions = BTreeMap::new();
    for f in files {
        if let Ok(bytes) = fs::read(&f) {
            if let Ok(v) = serde_json::from_slice::<Value>(&bytes) {
                scan(&v, &f, &mut additions);
            }
        }
    }
    let registry_path = "docs/curves/registry.json";
    let old = fs::read_to_string(registry_path).unwrap();
    let (start, end, raw) = segments(&old, "curves");
    let mut entries: Vec<(Value, String)> = raw
        .into_iter()
        .map(|s| (serde_json::from_str(&s).unwrap(), s))
        .collect();
    let mut browser_new = Vec::new();
    for (slug, (mut v, sources)) in additions {
        v["sources"] = json!(sources);
        if let Some((old, raw)) = entries.iter_mut().find(|(x, _)| x["slug"] == slug) {
            assert_eq!(old["icv1"], v["icv1"], "registry display slug collision");
            assert_eq!(
                old["model_json"], v["model_json"],
                "registry model identity collision"
            );
            let before = old.clone();
            let mut sources: Vec<String> = old["sources"]
                .as_array()
                .unwrap()
                .iter()
                .chain(v["sources"].as_array().unwrap())
                .map(|x| x.as_str().unwrap().to_string())
                .collect();
            sources.sort();
            sources.dedup();
            sources.truncate(6);
            old["sources"] = json!(sources);

            for s in v["sources"].as_array().unwrap() {
                let a = old["sources"].as_array_mut().unwrap();
                if a.len() < 6 && !a.contains(s) {
                    a.push(s.clone());
                }
            }
            for r in v["representations"].as_array().unwrap() {
                if !old["representations"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|x| x["curve_uid"] == r["curve_uid"])
                {
                    old["representations"]
                        .as_array_mut()
                        .unwrap()
                        .push(r.clone());
                }
            }
            if !old["representations"].as_array().unwrap().is_empty() {
                old.as_object_mut().unwrap().remove("ec1_unresolved");
            }
            if before != *old {
                *raw = render(old, 2);
            }
        } else {
            browser_new.push(v.clone());
            let raw = render(&v, 2);
            entries.push((v, raw));
        }
    }
    fn key(v: &Value) -> (u8, u64, BigUint, u64, String) {
        let fam = v["family"].as_str().unwrap();
        let rank = match fam {
            "koblitz" => 0,
            "subfield" => 1,
            "binary" => 2,
            "prime" => 3,
            _ => 4,
        };
        let p = v["params"]["p"]
            .as_str()
            .unwrap_or("0")
            .parse()
            .unwrap_or_default();
        (
            rank,
            v["params"]["n"]
                .as_u64()
                .or(v["params"]["m"].as_u64())
                .unwrap_or(0),
            p,
            if rank == 0 {
                v["params"]["a"].as_u64().unwrap_or(0)
            } else {
                0
            },
            v["slug"].as_str().unwrap().into(),
        )
    }
    entries.sort_by_key(|(v, _)| key(v));
    let text = format!(
        "{}\n{}\n {}",
        &old[..=start],
        entries
            .iter()
            .map(|(_, s)| format!("  {s}"))
            .collect::<Vec<_>>()
            .join(",\n"),
        &old[end..]
    );
    fs::write(registry_path, &text).unwrap();
    let ap = "src/cryptanalysis/curve_aliases.json";
    let mut aliases = read(ap);
    for (v, _) in &entries {
        for name in v["aliases"]
            .as_array()
            .unwrap()
            .iter()
            .chain(v["standard_names"].as_array().unwrap())
            .chain([&v["slug"], &v["icv1"]])
        {
            aliases["aliases"][curve_id::normalise(name.as_str().unwrap())] = v["slug"].clone();
        }
    }
    let mut buf = Vec::new();
    let fmt = serde_json::ser::PrettyFormatter::with_indent(b"");
    let mut ser = serde_json::Serializer::with_formatter(&mut buf, fmt);
    use serde::Serialize;
    aliases["aliases"].serialize(&mut ser).unwrap();
    let alias_text=format!("{{\n\"schema_version\": 1,\n\"generated_from\": \"docs/curves/registry.json\",\n\"generated_by\": \"scripts/build_curve_registry.py\",\n\"aliases\": {}\n}}\n",String::from_utf8(buf).unwrap());
    if fs::read(ap).ok().as_deref() != Some(alias_text.as_bytes()) {
        fs::write(ap, alias_text).unwrap();
    }
    write(
        root.join("prepared/registered-curves.json"),
        &json!({"new_models":browser_new.len(),"registry_sha256":sha256_hex(text.as_bytes()),"curves":entries.iter().filter(|(v,_)|v["sources"].as_array().unwrap().iter().any(|s|s.as_str().unwrap().starts_with(ROOT))).map(|(v,_)|v).collect::<Vec<_>>()}),
    );
    println!(
        "registered {} new models; browser indexing remains a separate check",
        browser_new.len()
    );
}
