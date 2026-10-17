//! Native, read-only summary of the frozen different-j full-image comparison.

use std::collections::BTreeMap;
use std::env;
use std::fs;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

fn median_u64(mut values: Vec<u64>) -> u64 {
    assert_eq!(values.len(), 3);
    values.sort_unstable();
    values[1]
}

fn median_f64(mut values: Vec<f64>) -> f64 {
    assert_eq!(values.len(), 3);
    values.sort_by(|a, b| a.total_cmp(b));
    values[1]
}

fn cpu(row: &Value) -> u64 {
    row["total_cpu_ns"].as_u64().expect("complete CPU")
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(
        args.len(),
        3,
        "usage: different_j_full_image_summary RUN REPLAY"
    );
    let raw = fs::read(&args[1]).expect("raw run");
    let replay_bytes = fs::read(&args[2]).expect("independent replay");
    let replay: Value = serde_json::from_slice(&replay_bytes).expect("replay JSON");
    assert_eq!(replay["checked"], 18);
    assert_eq!(replay["base_witness_checks"], 8);
    assert!(replay["relation_witness_checks"].as_u64().unwrap() > 0);
    assert_eq!(replay["records_sha256"], sha256_hex(&raw));
    let mut rows = BTreeMap::<(u64, String), Value>::new();
    for line in raw
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
    {
        let row: Value = serde_json::from_slice(line).expect("run JSON row");
        assert_eq!(row["verified"], true);
        assert_eq!(row["recovered"], 41);
        let key = (
            row["round"].as_u64().unwrap(),
            row["arm"].as_str().unwrap().to_string(),
        );
        assert!(rows.insert(key, row).is_none(), "duplicate arm-round row");
    }
    assert_eq!(rows.len(), 18);
    let arms = [
        "source",
        "isogeny-image",
        "destination-direct",
        "full-image",
        "rho-negation",
        "rho-negation-aa",
    ];
    let mut arm_summary = BTreeMap::new();
    for arm in arms {
        let values: Vec<&Value> = (0..3)
            .map(|round| &rows[&(round, arm.to_string())])
            .collect();
        let setup_values: Vec<u64> = values
            .iter()
            .filter_map(|row| row["setup_cpu_ns"].as_u64())
            .collect();
        arm_summary.insert(
            arm,
            json!({
                "median_total_cpu_ns": median_u64(values.iter().map(|row| cpu(row)).collect()),
                "median_setup_cpu_ns": if setup_values.len() == 3 { Some(median_u64(setup_values)) } else { None },
                "total_cpu_ns": values.iter().map(|row| cpu(row)).collect::<Vec<_>>(),
                "setup_cpu_ns": values.iter().map(|row| row["setup_cpu_ns"].as_u64()).collect::<Vec<_>>(),
                "collect_cpu_ns": values.iter().map(|row| row["collect_cpu_ns"].as_u64()).collect::<Vec<_>>(),
                "solve_cpu_ns": values.iter().map(|row| row["solve_cpu_ns"].as_u64()).collect::<Vec<_>>(),
                "verify_cpu_ns": values.iter().map(|row| row["verify_cpu_ns"].as_u64()).collect::<Vec<_>>(),
                "pair_keys": values.iter().map(|row| row["pair_keys"].as_u64()).collect::<Vec<_>>(),
                "rank": values.iter().map(|row| row["rank"].as_u64()).collect::<Vec<_>>(),
                "trials": values.iter().map(|row| row["trials"].as_u64()).collect::<Vec<_>>(),
                "relations": values.iter().map(|row| row["relations"].as_u64()).collect::<Vec<_>>(),
            }),
        );
    }
    let mut round_rows = Vec::new();
    let mut full_over_rho = Vec::new();
    let mut full_over_direct = Vec::new();
    let mut setup_fraction = Vec::new();
    for round in 0..3 {
        let full = &rows[&(round, "full-image".to_string())];
        let direct = &rows[&(round, "destination-direct".to_string())];
        let rho = &rows[&(round, "rho-negation".to_string())];
        assert_eq!(full["full_image_points"], 1525);
        assert_eq!(full["points_scanned"], 4575);
        assert_eq!(full["selection_candidates"], 21);
        assert_eq!(full["rank"], 9);
        assert_eq!(full["pair_keys"], direct["pair_keys"]);
        let full_cpu = cpu(full);
        let direct_cpu = cpu(direct);
        let rho_cpu = cpu(rho);
        let setup_cpu = full["setup_cpu_ns"].as_u64().unwrap();
        full_over_rho.push(full_cpu as f64 / rho_cpu as f64);
        full_over_direct.push(full_cpu as f64 / direct_cpu as f64);
        setup_fraction.push(setup_cpu as f64 / full_cpu as f64);
        round_rows.push(json!({
            "round":round,
            "full_total_cpu_ns":full_cpu,
            "full_setup_cpu_ns":setup_cpu,
            "direct_total_cpu_ns":direct_cpu,
            "rho_total_cpu_ns":rho_cpu,
            "full_over_rho_ratio":full_over_rho.last().unwrap(),
            "full_over_direct_ratio":full_over_direct.last().unwrap(),
            "full_setup_fraction":setup_fraction.last().unwrap(),
            "full_trials":full["trials"],
            "full_relations":full["relations"],
            "direct_trials":direct["trials"],
            "direct_relations":direct["relations"]
        }));
    }
    println!(
        "{}",
        json!({
            "schema":"different_j_full_image_summary/v1",
            "raw_sha256":sha256_hex(&raw),
            "replay_sha256":sha256_hex(&replay_bytes),
            "arm_rounds":rows.len(),
            "independent_scalar_replays":replay["checked"],
            "independent_base_witnesses":replay["base_witness_checks"],
            "independent_relation_witnesses":replay["relation_witness_checks"],
            "arms":arm_summary,
            "paired_rounds":round_rows,
            "median_full_over_rho":median_f64(full_over_rho),
            "median_full_over_direct":median_f64(full_over_direct),
            "median_full_setup_fraction":median_f64(setup_fraction)
        })
    );
}
