//! Reuse every higher-degree independent fixture; emit separately labeled controls.
use serde_json::Value;
use std::{collections::BTreeMap, fmt::Write, fs, path::Path};
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 2);
    let v: Value = serde_json::from_slice(&fs::read(&a[0]).unwrap()).unwrap();
    let base = Path::new(&a[1]);
    let mut jobs: Vec<_> = v["records"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|r| r["field_degree"] != 6 && r["mode"] == "population")
        .collect();
    jobs.sort_by_key(|r| {
        (
            r["field_degree"].as_u64().unwrap(),
            r["seed"].as_str().unwrap().parse::<u64>().unwrap(),
        )
    });
    let mut s = String::new();
    let mut fields = BTreeMap::new();
    for r in &jobs {
        let c = r["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["role"] == "source")
            .unwrap();
        let p = r["characteristic"].as_str().unwrap();
        let k = r["field_degree"].as_u64().unwrap();
        let m = serde_json::to_string(&r["modulus_low_coefficients"]).unwrap();
        writeln!(
            s,
            "construct_source({p},{k},{m},{},{},{},{});",
            serde_json::to_string(&c["original_root_coordinates"]).unwrap(),
            c["group_order"].as_str().unwrap(),
            r["seed"].as_str().unwrap(),
            serde_json::to_string(&c["icv1"]).unwrap()
        )
        .unwrap();
        fields.insert((k, p.to_string()), m);
    }
    assert_eq!(jobs.len(), 128);
    s.push_str("print(\"CONSTRUCTION_COMPLETE|sources=128\");\n");
    assert!(!base.join("general_sources.gp").exists());
    fs::write(base.join("general_sources.gp"), s).unwrap();
    let mut controls = String::new();
    for ((k, p), m) in fields {
        for i in 1..=4 {
            writeln!(
                controls,
                "general_control({p},{k},{m},{});",
                202610092000u64 + i
            )
            .unwrap();
        }
    }
    controls.push_str("print(\"ALGEBRA_COMPLETE|generalized_quadratic_controls=16\");\n");
    fs::write(base.join("general_control_inputs.gp"), controls).unwrap();
    println!("128 original sources and 16 separate higher-degree quadratic controls emitted");
}
