//! Emit frozen ecbench specs from research/drastic_prime_conductor_gap_20261012/pairs/pairs.json.
//! One spec per pair: crater first, floor second. rho arms on every pair;
//! the ic-mitm arm (the prior round's exact configuration) only on the
//! smallest-ell pairs where the 3600 s budget can hold.
use serde_json::{json, Value};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pairs_path = args
        .iter()
        .position(|s| s == "--pairs")
        .and_then(|i| args.get(i + 1))
        .expect("--pairs <path>")
        .clone();
    let out_dir = args
        .iter()
        .position(|s| s == "--out")
        .and_then(|i| args.get(i + 1))
        .expect("--out <dir>")
        .clone();
    std::fs::create_dir_all(&out_dir).unwrap();

    let doc: Value =
        serde_json::from_slice(&std::fs::read(&pairs_path).unwrap()).expect("pairs json");
    let mut written = Vec::new();
    for (_i, pair) in doc["pairs"].as_array().unwrap().iter().enumerate() {
        let d = pair["d"].as_u64().unwrap();
        let ell = pair["ell"].as_u64().unwrap();
        let label = format!("d{}-ell{}", d, ell);
        let q_bits = 64 - pair["p"].as_u64().unwrap().leading_zeros();
        let small = ell < 600_000; // ic-mitm only on the smallest-ell pairs
        let arms = if small {
            json!([
                {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
                {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}},
                {"name": "ic-mitm", "role": "candidate", "method": {
                    "id": "ic.pipeline",
                    "params": {
                        "factor_base": "prime-abscissa:size=2048",
                        "oracle": "mitm:m=2,negation_folded=1",
                        "linalg": "structured-gauss",
                        "targets": "walk",
                        "max_trials": "100000000",
                        "solver_budget_seconds": "0"
                    }
                }}
            ])
        } else {
            json!([
                {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
                {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}}
            ])
        };
        let spec = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic prime conductor gap: D=-{}, ell={}, {}-bit field", d, ell, q_bits),
            "description": "Fresh certified CM pair from the drastic-gap search; independent one-target workloads on each curve, target law hashes the curve model. No isogeny walk; no cross-curve transport. Counts are primary.",
            "workloads": {
                "curves": [pair["crater"].clone(), pair["floor"].clone()],
                "targets_per_curve": 8,
                "target_seed": 2026101201
            },
            "arms": arms,
            "measurement": {
                "rounds": 1,
                "warmup": 1,
                "order": "alternate",
                "seed": 20261012,
                "isolation_required": "L2",
                "timeout_seconds": if small { 3600 } else { 600 }
            }
        });
        let path = format!("{}/{}.json", out_dir, label);
        std::fs::write(&path, serde_json::to_vec_pretty(&spec).unwrap()).unwrap();
        written.push(path);
    }
    for w in &written {
        println!("wrote {}", w);
    }
}
