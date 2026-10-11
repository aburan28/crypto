//! Preserve original independent source fixtures when changing the search policy.
use serde_json::Value;
use std::{fmt::Write, fs};
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 2, "SAVED_RECORDS OUTPUT_GP");
    let v: Value = serde_json::from_slice(&fs::read(&a[0]).unwrap()).unwrap();
    let mut records: Vec<_> = v["records"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|r| r["mode"] == "population" && r["field_degree"] == 6)
        .collect();
    records.sort_by_key(|r| {
        (
            r["characteristic"]
                .as_str()
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            r["seed"].as_str().unwrap().parse::<u64>().unwrap(),
        )
    });
    let mut s = String::new();
    for r in &records {
        let src = r["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["role"] == "source")
            .unwrap();
        writeln!(
            s,
            "construct_source({},{},{},{},{},{});",
            r["characteristic"].as_str().unwrap(),
            serde_json::to_string(&r["modulus_low_coefficients"]).unwrap(),
            serde_json::to_string(&src["original_root_coordinates"]).unwrap(),
            src["group_order"].as_str().unwrap(),
            r["seed"].as_str().unwrap(),
            serde_json::to_string(&src["icv1"]).unwrap()
        )
        .unwrap();
    }
    writeln!(
        s,
        "print(\"CONSTRUCTION_COMPLETE|sources={}\");",
        records.len()
    )
    .unwrap();
    assert!(
        !std::path::Path::new(&a[1]).exists(),
        "preserve previous input freeze"
    );
    fs::write(&a[1], s).unwrap();
    println!("emitted {} degree-6 independent sources", records.len());
}
