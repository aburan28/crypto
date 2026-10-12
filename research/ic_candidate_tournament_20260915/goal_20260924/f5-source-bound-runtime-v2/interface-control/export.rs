//! Deserialization control for the exact pinned native Job and Config types.
#[path = "ic_tournament_worker.rs"]
mod pinned_worker;
use serde_json::{json, Value};
use std::io::Read;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let cases: Vec<Value> = serde_json::from_str(&input).unwrap();
    let results: Vec<Value> = cases.iter().map(|case| {
        match pinned_worker::inspect_job_schema(&case["input"].to_string()) {
            Ok(observed) => json!({"label": case["label"], "status": "PARSED", "observed": observed}),
            Err(error) => json!({"label": case["label"], "status": "REJECTED", "error": error})
        }
    }).collect();
    println!("{}", json!({"schema_version": 1, "status": "NATIVE_JOB_SCHEMA_CONTROL", "results": results}));
}
