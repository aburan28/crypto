//! Known-answer, one-map selection from a previously frozen verified record.
use serde_json::{json, Value};
use std::{fs, path::Path};
#[path = "../../src/hash/sha256.rs"]
mod sha;
fn hash(b: &[u8]) -> String {
    sha::sha256(b).iter().map(|x| format!("{x:02x}")).collect()
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 2);
    let bytes = fs::read(&a[0]).unwrap();
    let mut doc: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(doc["status"], "PASS");
    let first = doc["result"]["isogenies"][0].clone();
    doc["result"]["isogenies"] = json!([first]);
    doc["result"]["count"] = json!(1);
    doc["scope"] = json!("one_frobenius_eigenline");
    doc["degree_coverage"] = json!("PARTIAL");
    doc["expected_total_eigenlines"] = json!(2);
    doc["constructed_eigenlines"] = json!(1);
    doc["known_answer_control"] = json!({"source":a[0],"sha256":hash(&bytes),"selection":"first previously verified map; no new construction"});
    let out = Path::new(&a[1]);
    fs::create_dir_all(out.join("p192")).unwrap();
    let output = serde_json::to_vec_pretty(&doc).unwrap();
    fs::write(out.join("p192/ell-199.json"), &output).unwrap();
    fs::write(out.join("p192/ell-199.stderr.txt"), b"").unwrap();
    let receipt = json!({"command":["known-answer-map-selection",a[0]],"status":"PASS","exit_code":0,"stdout":"ell-199.json","stdout_sha256":hash(&output),"stderr_sha256":hash(b""),"map_count":1,"known_answer_control":true});
    fs::write(
        out.join("p192/ell-199.receipt.json"),
        serde_json::to_vec_pretty(&receipt).unwrap(),
    )
    .unwrap();
    fs::write(out.join("search.json"),serde_json::to_vec_pretty(&json!({"attempts":[{"curve":"p192","ell":199,"receipt":receipt}],"known_answer_control":true})).unwrap()).unwrap();
}
