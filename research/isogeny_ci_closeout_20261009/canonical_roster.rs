//! Restore canonical roster key order without recalculating or changing evidence.
use serde::{
    de::{MapAccess, Visitor},
    Deserialize, Deserializer,
};
use serde_json::{value::RawValue, Value};
use std::{fmt, fs};
#[path = "../../src/hash/sha256.rs"]
#[allow(dead_code)]
mod sha;
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
    fn raw(&self, key: &str) -> &str {
        self.0.iter().find(|(k, _)| k == key).unwrap().1.get()
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        (2..=3).contains(&args.len()),
        "canonical-roster FILE [--check]"
    );
    let path = &args[1];
    let before = fs::read_to_string(path).unwrap();
    let reference: Value = serde_json::from_str(&before).unwrap();
    let mut doc: Object = serde_json::from_str(&before).unwrap();
    if args.get(2).is_some_and(|arg| arg == "--browser-sources") {
        let mut sources: Object = serde_json::from_str(doc.raw("sources")).unwrap();
        let mut expected = reference;
        for source in [
            "docs/curves/registry.json",
            "docs/curves/covers.json",
            "docs/curves/cover-links.yaml",
            "docs/ic/leaderboard.json",
        ] {
            let checksum: String = sha::sha256(&fs::read(source).unwrap())
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            sources.0.iter_mut().find(|(k, _)| k == source).unwrap().1 =
                RawValue::from_string(serde_json::to_string(&checksum).unwrap()).unwrap();
            expected["sources"][source] = Value::String(checksum);
        }
        let value = format!(
            "{{\n  {}\n }}",
            sources
                .0
                .iter()
                .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
                .collect::<Vec<_>>()
                .join(",\n  ")
        );
        doc.0.iter_mut().find(|(k, _)| k == "sources").unwrap().1 =
            RawValue::from_string(value).unwrap();
        let after = format!(
            "{{\n {}\n}}\n",
            doc.0
                .iter()
                .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
                .collect::<Vec<_>>()
                .join(",\n ")
        );
        assert_eq!(
            serde_json::from_str::<Value>(&after).unwrap(),
            expected,
            "browser content changed beyond its source hashes"
        );
        fs::write(path, after).unwrap();
        println!("Browser source hashes refreshed; every other value preserved");
        return;
    }
    let board = doc.raw("board").to_owned();
    let rows: Vec<Box<RawValue>> = serde_json::from_str(doc.raw("roster")).unwrap();
    let order = [
        "slug",
        "family",
        "field",
        "coefficients",
        "order_bits",
        "standard",
        "ec1",
        "ec1_unresolved",
        "legacy",
        "on_board",
        "curves_yaml_key",
        "factor_base_link_status",
    ];
    let mut formatted = vec![];
    for row in rows {
        let object: Object = serde_json::from_str(row.get()).unwrap();
        assert_eq!(
            object.0.len(),
            order.len(),
            "unknown roster field must be reviewed"
        );
        let fields: Vec<_> = order
            .iter()
            .map(|key| {
                format!(
                    "{}: {}",
                    serde_json::to_string(key).unwrap(),
                    object.raw(key)
                )
            })
            .collect();
        formatted.push(format!("{{\n   {}\n  }}", fields.join(",\n   ")));
    }
    let roster = format!("[\n  {}\n ]", formatted.join(",\n  "));
    doc.0.iter_mut().find(|(k, _)| k == "roster").unwrap().1 =
        RawValue::from_string(roster).unwrap();
    assert_eq!(doc.raw("board"), board, "measurement bytes changed");
    let after = format!(
        "{{\n {}\n}}\n",
        doc.0
            .iter()
            .map(|(k, v)| format!("{}: {}", serde_json::to_string(k).unwrap(), v.get()))
            .collect::<Vec<_>>()
            .join(",\n ")
    );
    assert_eq!(
        serde_json::from_str::<Value>(&after).unwrap(),
        reference,
        "a value changed"
    );
    if args.len() == 3 {
        assert_eq!(args[2], "--check");
        assert_eq!(after, before, "roster serialization is stale");
        println!(
            "PASS: canonical roster serialization; every value and measurement byte preserved"
        );
    } else {
        fs::write(path, after).unwrap();
        println!("Canonical roster serialization updated; all values preserved");
    }
}
