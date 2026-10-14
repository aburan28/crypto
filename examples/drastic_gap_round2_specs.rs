//! Emit frozen ecbench specs for the round-2 drastic-gap study. Two specs
//! per D = -3 pair (shared target seed, so the crater's targets are
//! identical): a rho spec on both curves, and a crater-only rho.cubic
//! spec. Plus the round-1 lead replication (d43-ell1016839, 32 fresh
//! targets). 32 targets per curve; the cubic arm refuses the floor by
//! construction, so it runs in its own spec.
use serde_json::{json, Value};
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pairs_path = args
        .iter()
        .position(|s| s == "--pairs")
        .and_then(|i| args.get(i + 1))
        .expect("--pairs <path>")
        .clone();
    let round1_path = args
        .iter()
        .position(|s| s == "--round1")
        .and_then(|i| args.get(i + 1))
        .expect("--round1 <path>")
        .clone();
    let out_dir = args
        .iter()
        .position(|s| s == "--out")
        .and_then(|i| args.get(i + 1))
        .expect("--out <dir>")
        .clone();
    std::fs::create_dir_all(&out_dir).unwrap();

    let rho_arms = json!([
        {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
        {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}}
    ]);
    let cubic_arms = json!([
        {"name": "rho-cubic", "role": "candidate", "method": {"id": "rho.cubic"}}
    ]);
    let measurement = |timeout: u64| json!({
        "rounds": 1, "warmup": 1, "order": "alternate",
        "seed": 20261012, "isolation_required": "L2",
        "timeout_seconds": timeout
    });

    let doc: Value =
        serde_json::from_slice(&std::fs::read(&pairs_path).unwrap()).expect("pairs json");
    let mut written = Vec::new();
    for pair in doc["pairs"].as_array().unwrap() {
        let ell = pair["ell"].as_u64().unwrap();
        let label = format!("d3-ell{}", ell);
        let rho = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 2, rho: D=-3, ell={}", ell),
            "description": "Certified j = 0 crater and j != 0 floor of a D = -3 pair; independent one-target workloads. Counts are primary.",
            "workloads": {"curves": [pair["crater"].clone(), pair["floor"].clone()],
                          "targets_per_curve": 32, "target_seed": 2026101202},
            "arms": rho_arms.clone(),
            "measurement": measurement(600),
        });
        let cubic = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 2, cubic crater: D=-3, ell={}", ell),
            "description": "The same crater as the rho spec (same target seed, so identical targets), walked by rho.cubic over its six automorphism classes.",
            "workloads": {"curves": [pair["crater"].clone()],
                          "targets_per_curve": 32, "target_seed": 2026101202},
            "arms": cubic_arms.clone(),
            "measurement": measurement(600),
        });
        for (suffix, spec) in [("rho", rho), ("cubic", cubic)] {
            let path = format!("{}/{}-{}.json", out_dir, label, suffix);
            std::fs::write(&path, serde_json::to_vec_pretty(&spec).unwrap()).unwrap();
            written.push(path);
        }
    }

    // round-1 lead replication: the d43-ell1016839 pair, 32 fresh targets
    let r1: Value =
        serde_json::from_slice(&std::fs::read(&round1_path).unwrap()).expect("round-1 pairs");
    let lead = r1["pairs"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["d"] == 43 && p["ell"] == 1016839)
        .expect("the round-1 lead pair");
    let repl = json!({
        "schema": "ecbench.spec/v1",
        "label": "Drastic gap round 2 replication: d43-ell1016839, 32 fresh targets",
        "description": "Round 1 read a 95%-level 0.5286 rho floor/crater quotient on this pair with 8 workloads; 32 fresh independent targets re-measure it. No isogeny, no transport. Counts are primary.",
        "workloads": {"curves": [lead["crater"].clone(), lead["floor"].clone()],
                      "targets_per_curve": 32, "target_seed": 2026101203},
        "arms": rho_arms.clone(),
        "measurement": measurement(900),
    });
    let path = format!("{}/replicate-d43-ell1016839.json", out_dir);
    std::fs::write(&path, serde_json::to_vec_pretty(&repl).unwrap()).unwrap();
    written.push(path);

    for w in &written {
        println!("wrote {}", w);
    }
}
