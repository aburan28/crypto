//! Validate the finite coverage and stream roster, then summarize recorded counts.
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
};
fn main() {
    let path = std::env::args().nth(1).expect("IDENTIFIED.jsonl");
    let text = fs::read_to_string(path).unwrap();
    let mut curves = BTreeSet::new();
    let mut uids = BTreeSet::new();
    let mut degrees = BTreeMap::<u64, u64>::new();
    let mut seen = BTreeSet::new();
    let (mut checks, mut runs, mut trials, mut relations, mut solved, mut capped) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    for line in text.lines() {
        let v: Value = serde_json::from_str(line).unwrap();
        match v["kind"].as_str().unwrap() {
            "manifest" => {
                assert_eq!(v["vertical_Fp_edges_exist"], false);
                assert_eq!(v["max_targets_per_run"], 100000);
                assert_eq!(v["runs_per_curve"], 12);
            }
            "curve" => {
                assert!(curves.insert(v["edge_id"].as_str().unwrap().to_string()));
                assert!(uids.insert(v["measured_curve_uid"].as_str().unwrap().to_string()));
                assert!(v["measured_curve_uid"]
                    .as_str()
                    .unwrap()
                    .starts_with("urn:ec-record:1:sha256:"));
            }
            "transport_check" => {
                assert_eq!(v["status"], "PASS");
                assert_eq!(v["scalar_checks"], 64);
                assert_eq!(v["transported_pair_sums_checked"], 36);
                checks += 1;
                *degrees.entry(v["degree"].as_u64().unwrap()).or_default() += 1;
            }
            "run" => {
                let id = v["edge_id"].as_str().unwrap().to_string();
                let run = v["run"].as_u64().unwrap();
                assert!(run < 12);
                assert_eq!(v["seed"].as_u64().unwrap(), 2026101600 + run);
                assert!(seen.insert((id, run)));
                runs += 1;
                trials += v["targets"].as_u64().unwrap();
                relations += v["relations"].as_u64().unwrap();
                assert_eq!(v["rejected"], 0);
                if v["status"] == "EXHAUSTED_TRIAL_CAP" {
                    assert_eq!(v["targets"], 100000);
                    assert_eq!(v["verified"], false);
                    capped += 1;
                }
                if v["status"] == "VERIFIED_SOLVE" {
                    assert_eq!(v["verified"], true);
                    solved += 1;
                }
            }
            _ => panic!("unknown record type"),
        }
    }
    assert_eq!(
        degrees,
        BTreeMap::from([
            (3, 1),
            (7, 2),
            (11, 2),
            (13, 2),
            (19, 2),
            (23, 2),
            (53, 2),
            (61, 2)
        ])
    );
    assert_eq!(checks, 15);
    assert_eq!(curves.len(), 16);
    assert_eq!(runs, 192);
    for id in &curves {
        for run in 0..12 {
            assert!(seen.contains(&(id.clone(), run)));
        }
    }
    println!("{}",serde_json::to_string_pretty(&json!({"schema":"prime61.isogeny-summary/v1","outgoing_horizontal_kernel_edges":checks,
        "source_and_destination_workloads":curves.len(),"prime_degree_bound":61,"degree_histogram":degrees,
        "runs_per_workload":12,"runs":runs,"target_trials":trials,"relations":relations,"verified_solves":solved,"capped_runs":capped,
        "transport_scalar_checks":checks*64,"transported_pair_sum_checks":checks*36,
        "vertical_Fp_edges_exist":false,"all_horizontal_class_vertices_exhausted":false,
        "full_operation_cost":null,"speedup":null})).unwrap());
}
