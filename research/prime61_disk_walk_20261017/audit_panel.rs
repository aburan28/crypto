//! Independent roster/counter comparison against the frozen one-hop JSONL.
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader},
};
fn rows(path: &str) -> impl Iterator<Item = Value> {
    BufReader::new(File::open(path).unwrap())
        .lines()
        .map(|s| serde_json::from_str(&s.unwrap()).unwrap())
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    assert_eq!(a.len(), 3, "OLD_IDENTIFIED.jsonl NEW_CHECKPOINT.jsonl");
    let mut old_ids = BTreeMap::new();
    let mut old_points = BTreeMap::new();
    let mut wanted = BTreeMap::new();
    for v in rows(&a[1]) {
        match v["kind"].as_str().unwrap() {
            "curve" => {
                let uid = v["measured_curve_uid"].as_str().unwrap().to_string();
                old_ids.insert(v["edge_id"].as_str().unwrap().to_string(), uid.clone());
                old_points.insert(
                    uid,
                    (v["transported_G"].clone(), v["transported_Q"].clone()),
                );
            }
            "run" => {
                let uid = old_ids[v["edge_id"].as_str().unwrap()].clone();
                let run = v["run"].as_u64().unwrap();
                wanted.insert((uid,run),json!({"seed":v["seed"],"targets":v["targets"],"relations":v["relations"],"scalar":v["scalar"],"status":v["status"]}));
            }
            _ => (),
        }
    }
    assert_eq!(wanted.len(), 192);
    let mut new_ids = BTreeMap::new();
    let mut matches = 0;
    for v in rows(&a[2]) {
        match v["kind"].as_str().unwrap() {
            "curve" if v["id"].as_i64().unwrap() <= 16 => {
                let p = &v["payload"];
                let uid = p["uid"].as_str().unwrap().to_string();
                let point = |k: &str| {
                    Value::Array(
                        p[k].as_array()
                            .unwrap()
                            .iter()
                            .map(|x| Value::String(x.as_u64().unwrap().to_string()))
                            .collect(),
                    )
                };
                assert_eq!(old_points[&uid], (point("g"), point("q")));
                new_ids.insert(v["id"].as_i64().unwrap(), uid);
            }
            "run" if v["curve"].as_i64().unwrap() <= 16 => {
                let p = &v["payload"];
                let uid = new_ids[&v["curve"].as_i64().unwrap()].clone();
                let run = v["run"].as_u64().unwrap();
                let got = json!({"seed":p["seed"],"targets":p["targets"],"relations":p["relations"].as_array().unwrap().len(),"scalar":p["scalar"].as_u64().map(|x|x.to_string()),"status":p["status"]});
                assert_eq!(wanted.remove(&(uid, run)).unwrap(), got);
                matches += 1;
            }
            _ => (),
        }
    }
    assert!(wanted.is_empty());
    assert_eq!(new_ids.len(), 16);
    assert_eq!(matches, 192);
    println!(
        "{}",
        json!({"schema":"prime61.disk-baseline-audit/v1","status":"PASS","exact_transported_workloads_matched":16,"paired_run_rows_matched":matches,"source":"frozen degree-61 panel","full_cost_comparison":null})
    );
}
