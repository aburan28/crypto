//! Native, incremental catalogue joins. Existing measurement records retain their bytes.
use num_bigint::BigUint;
use serde::{
    de::{MapAccess, Visitor},
    Deserialize, Deserializer,
};
use serde_json::{json, value::RawValue, Value};
use std::{collections::BTreeMap, fmt, fs};
#[path = "../../src/hash/sha256.rs"]
#[allow(dead_code)]
mod sha;
fn digest(b: &[u8]) -> String {
    sha::sha256(b).iter().map(|b| format!("{b:02x}")).collect()
}
fn load(path: &str) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn integer(v: &Value) -> BigUint {
    BigUint::parse_bytes(v.as_str().unwrap().as_bytes(), 10).unwrap()
}
fn rendered(v: &Value, indent: usize) -> String {
    let mut bytes = Vec::new();
    let mut serializer = serde_json::Serializer::with_formatter(
        &mut bytes,
        serde_json::ser::PrettyFormatter::with_indent(b" "),
    );
    serde::Serialize::serialize(v, &mut serializer).unwrap();
    String::from_utf8(bytes)
        .unwrap()
        .lines()
        .enumerate()
        .map(|(i, l)| {
            if i == 0 {
                l.into()
            } else {
                format!("{}{l}", " ".repeat(indent))
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}
struct Object(Vec<(String, Box<RawValue>)>);
impl<'de> Deserialize<'de> for Object {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct Ordered;
        impl<'de> Visitor<'de> for Ordered {
            type Value = Object;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("JSON object")
            }
            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Object, A::Error> {
                let mut fields = vec![];
                while let Some(field) = map.next_entry()? {
                    fields.push(field);
                }
                Ok(Object(fields))
            }
        }
        d.deserialize_map(Ordered)
    }
}
impl Object {
    fn read(path: &str) -> Self {
        serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
    }
    fn raw(&self, key: &str) -> &str {
        self.0.iter().find(|(k, _)| k == key).unwrap().1.get()
    }
    fn value(&self, key: &str) -> Value {
        serde_json::from_str(self.raw(key)).unwrap()
    }
    fn set_raw(&mut self, key: &str, s: String) {
        self.0.iter_mut().find(|(k, _)| k == key).unwrap().1 = RawValue::from_string(s).unwrap();
    }
    fn set(&mut self, key: &str, v: &Value) {
        self.set_raw(key, rendered(v, 1));
    }
    fn write(&self, path: &str) {
        let text = format!(
            "{{\n {}\n}}\n",
            self.0
                .iter()
                .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
                .collect::<Vec<_>>()
                .join(",\n ")
        );
        fs::write(path, text).unwrap();
    }
}
fn raw_rows(raw: &str) -> BTreeMap<String, Box<RawValue>> {
    let rows: Vec<Box<RawValue>> = serde_json::from_str(raw).unwrap();
    rows.into_iter()
        .map(|r| {
            let v: Value = serde_json::from_str(r.get()).unwrap();
            (v["slug"].as_str().unwrap().into(), r)
        })
        .collect()
}
fn array(rows: Vec<String>) -> String {
    format!("[\n  {}\n ]", rows.join(",\n  "))
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#x27;")
}
fn main() {
    assert_eq!(
        std::env::args().len(),
        1,
        "run catalogue-views from the repository root"
    );
    let registry = load("docs/curves/registry.json");
    let covers = load("docs/curves/covers.json");
    let graph = load("docs/curves/cover-links.yaml");
    let registry_hash = digest(&fs::read("docs/curves/registry.json").unwrap());
    assert_eq!(covers["registry_sha256"], registry_hash);
    assert_eq!(graph["registry_sha256"], registry_hash);
    assert_eq!(covers["summary"]["invalid_input"], 0);
    let curves = registry["curves"].as_array().unwrap();
    let mut leaderboard = Object::read("docs/ic/leaderboard.json");
    let board_bytes = leaderboard.raw("board").to_owned();
    let roster_old = raw_rows(leaderboard.raw("roster"));
    let old_count = roster_old.len();
    let old_hash = leaderboard.value("sources")["registry"]["sha256"]
        .as_str()
        .unwrap()
        .to_owned();
    let mut roster = vec![];
    for c in curves {
        let slug = c["slug"].as_str().unwrap();
        if let Some(old) = roster_old.get(slug) {
            roster.push(old.get().into());
            continue;
        }
        assert_eq!(
            c["family"], "prime",
            "new model outside this study's prime-field scope"
        );
        let reps = c["representations"].as_array().unwrap();
        let row = json!({"slug":slug,"family":"prime","field":format!("GF(p), {}-bit p",integer(&c["params"]["p"]).bits()),
            "coefficients":"standard coefficients","order_bits":integer(&c["order"]).bits(),"standard":c["standard_names"],
            "ec1":reps.iter().map(|r|r["ec1"].clone()).collect::<Vec<_>>(),"ec1_unresolved":null,"legacy":c["aliases"],
            "on_board":false,"curves_yaml_key":null,"factor_base_link_status":"no_curves_yaml_record"});
        roster.push(rendered(&row, 2));
    }
    leaderboard.set_raw("roster", array(roster));
    let raw = leaderboard
        .raw("sources")
        .replace(&old_hash, &registry_hash);
    leaderboard.set_raw("sources", raw);
    assert_eq!(
        leaderboard.raw("board"),
        board_bytes,
        "measurement bytes changed"
    );
    let roster = leaderboard.value("roster");
    let rows = roster.as_array().unwrap();
    let prime: Vec<_> = rows.iter().filter(|r| r["family"] == "prime").collect();
    let mut html = fs::read_to_string("docs/ic-leaderboard.html").unwrap();
    let start = html
        .find("<details ><summary>prime ·")
        .expect("prime roster section");
    let end = start + html[start..].find("</details>").unwrap() + "</details>".len();
    let mut section=format!("<details ><summary>prime · {} curves · {} on the board</summary><div class=\"scroll\"><table><thead><tr><th>slug</th><th>field</th><th>coefficients</th><th class=\"n\">#E bits</th><th>EC1</th><th>also called</th><th>curves.yaml</th><th>board</th></tr></thead><tbody>\n",prime.len(),prime.iter().filter(|r|r["on_board"]==true).count());
    for c in prime {
        let ids = c["ec1"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>()
            .join(", ");
        let ec1 = if ids.is_empty() {
            format!(
                "<code></code><span class=\"muted\">{}</span>",
                escape(c["ec1_unresolved"].as_str().unwrap_or("None"))
            )
        } else {
            format!("<code>{}</code>", escape(&ids))
        };
        let aliases = c["standard"]
            .as_array()
            .unwrap()
            .iter()
            .chain(c["legacy"].as_array().unwrap())
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>()
            .join(", ");
        section.push_str(&format!("<tr><td><span class=\"slug\">{}</span></td><td>{}</td><td>{}</td><td class=\"n\">{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",escape(c["slug"].as_str().unwrap()),escape(c["field"].as_str().unwrap()),escape(c["coefficients"].as_str().unwrap()),c["order_bits"],ec1,escape(&aliases),escape(c["curves_yaml_key"].as_str().unwrap_or("—")),if c["on_board"]==true {"✓"} else {""}));
    }
    section.push_str("</tbody></table></div></details>");
    html.replace_range(start..end, &section);
    html = html.replace(
        &format!("of {old_count} named in the repository"),
        &format!("of {} named in the repository", curves.len()),
    );
    html = html.replace(
        &format!(
            "<code>docs/curves/registry.json</code> ({})",
            &old_hash[..12]
        ),
        &format!(
            "<code>docs/curves/registry.json</code> ({})",
            &registry_hash[..12]
        ),
    );
    let md = fs::read_to_string("docs/ic/LEADERBOARD.md")
        .unwrap()
        .replace(
            &format!(
                "`docs/curves/registry.json` — sha256 `{}…`",
                &old_hash[..16]
            ),
            &format!(
                "`docs/curves/registry.json` — sha256 `{}…`",
                &registry_hash[..16]
            ),
        );
    leaderboard.write("docs/ic/leaderboard.json");
    fs::write("docs/ic-leaderboard.html", html).unwrap();
    fs::write("docs/ic/LEADERBOARD.md", md).unwrap();
    let mut aliases = Object::read("src/cryptanalysis/curve_aliases.json");
    let mut names = serde_json::Map::new();
    for c in curves {
        for value in c["aliases"]
            .as_array()
            .unwrap()
            .iter()
            .chain(c["standard_names"].as_array().unwrap())
            .chain([&c["slug"], &c["icv1"]])
        {
            let key = value
                .as_str()
                .unwrap()
                .trim()
                .trim_matches('`')
                .chars()
                .filter(|ch| !ch.is_whitespace() && !matches!(ch, '{' | '}' | '_'))
                .map(|ch| match ch {
                    '₀'..='₉' => char::from_u32('0' as u32 + ch as u32 - '₀' as u32).unwrap(),
                    _ => ch,
                })
                .collect::<String>()
                .to_lowercase();
            names.insert(key, c["slug"].clone());
        }
    }
    aliases.set_raw(
        "aliases",
        format!(
            "{{\n{}\n}}",
            names
                .iter()
                .map(|(k, v)| format!(
                    "{}: {}",
                    serde_json::to_string(k).unwrap(),
                    serde_json::to_string(v).unwrap()
                ))
                .collect::<Vec<_>>()
                .join(",\n")
        ),
    );
    // The historical alias map uses indentation zero.
    let alias_text = format!(
        "{{\n{}\n}}\n",
        aliases
            .0
            .iter()
            .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
            .collect::<Vec<_>>()
            .join(",\n")
    );
    fs::write("src/cryptanalysis/curve_aliases.json", alias_text).unwrap();
    let cover_by_slug: BTreeMap<_, _> = covers["curves"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| (r["slug"].as_str().unwrap(), r))
        .collect();
    let nodes: BTreeMap<_, _> = graph["curves"]
        .as_array()
        .unwrap()
        .iter()
        .map(|r| (r["icv1_identity"]["slug"].as_str().unwrap(), r))
        .collect();
    assert_eq!(cover_by_slug.len(), curves.len());
    assert_eq!(nodes.len(), curves.len());
    let mut browser = Object::read("docs/browser/data.json");
    let existing = raw_rows(browser.raw("curves"));
    let mut joined = vec![];
    for c in curves {
        let slug = c["slug"].as_str().unwrap();
        let mut row: Value = if let Some(old) = existing.get(slug) {
            serde_json::from_str(old.get()).unwrap()
        } else {
            let p = integer(&c["params"]["p"]);
            let rep = &c["representations"][0];
            let subgroup = &rep["curve"];
            json!({"slug":slug,"icv1":c["icv1"],"family":"prime","standard_names":c["standard_names"],"field":{"characteristic":"p","p":p.to_string(),"bits":p.bits(),"degree":1,"modulus":null},
                "a":c["params"]["a"],"b":c["params"]["b"],"construction":null,"trace":c["trace"].to_string(),"order":c["order"],"order_bits":integer(&c["order"]).bits(),
                "r":subgroup["subgroup_order"],"r_bits":integer(&subgroup["subgroup_order"]).bits(),"cofactor":subgroup["cofactor"],"generator":subgroup["generator"],"target_group":subgroup["target_group"],
                "j":c["j"],"endomorphism_discriminant":null,"model_json":c["model_json"],"ec1":[rep["ec1"]],"curve_uid":[rep["curve_uid"]],"ec1_unresolved":null,"retired_names":c["aliases"],"sources":c["sources"],
                "on_leaderboard":false,"leaderboard":[],"ecbench":[],"factor_bases":[],"tournament_cells":[]})
        };
        assert_eq!(
            cover_by_slug[slug]["model_sha256"],
            digest(c["model_json"].as_str().unwrap().as_bytes())
        );
        assert_eq!(nodes[slug]["model_json"], c["model_json"]);
        row["hyperelliptic_cover"] = cover_by_slug[slug].clone();
        row["standards_provenance"] = nodes[slug]["standards_provenance"]
            .as_array()
            .map_or(json!([]), |a| json!(a));
        row["cover_links"] = json!(nodes[slug]["cover_links"]
            .as_array()
            .unwrap()
            .iter()
            .map(|link| {
                let mut l = link.clone();
                l["cover_id"] =
                    graph["covers"][link["cover_uid"].as_str().unwrap()]["cover_id"].clone();
                l["map_id"] = graph["maps"][link["map_uid"].as_str().unwrap()]["map_id"].clone();
                l
            })
            .collect::<Vec<_>>());
        if let Some(old) = existing.get(slug) {
            if serde_json::from_str::<Value>(old.get()).unwrap() == row {
                joined.push(old.get().into());
                continue;
            }
        }
        joined.push(rendered(&row, 2));
    }
    browser.set_raw("curves", array(joined));
    let mut counts: Object = serde_json::from_str(browser.raw("counts")).unwrap();
    counts.set_raw("curves", curves.len().to_string());
    browser.set_raw(
        "counts",
        format!(
            "{{\n  {}\n }}",
            counts
                .0
                .iter()
                .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
                .collect::<Vec<_>>()
                .join(",\n  ")
        ),
    );
    let mut sources: Object = serde_json::from_str(browser.raw("sources")).unwrap();
    for path in [
        "docs/curves/registry.json",
        "docs/curves/covers.json",
        "docs/curves/cover-links.yaml",
        "docs/ic/leaderboard.json",
    ] {
        sources.set_raw(
            path,
            serde_json::to_string(&digest(&fs::read(path).unwrap())).unwrap(),
        );
    }
    browser.set_raw(
        "sources",
        format!(
            "{{\n  {}\n }}",
            sources
                .0
                .iter()
                .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
                .collect::<Vec<_>>()
                .join(",\n  ")
        ),
    );
    browser.write("docs/browser/data.json");
    println!("Catalogue joins refreshed: {} models; {} new roster records; existing board bytes preserved.",curves.len(),curves.len()-old_count);
}
