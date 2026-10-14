//! Emit frozen ecbench specs for round 3. All specs use the paired
//! planted_shared target law (same scalars on both curves of a pair).
//! Per pair: a rho spec on both curves and an automorphism spec on the
//! crater only (same target seed, so the crater solves the same ks).
//! Plus the D = -19 ell ~ 2^23 gap-sweep pair (rho only, 32 targets).
use serde_json::{json, Value};
use std::path::Path;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pairs_path = arg(&args, "--pairs");
    let out_dir = arg(&args, "--out");
    std::fs::create_dir_all(&out_dir).unwrap();

    let doc: Value =
        serde_json::from_slice(&std::fs::read(&pairs_path).unwrap()).expect("pairs json");
    // round-2 D = -3 pairs re-measured under the paired law with 64 targets
    let round2_path = arg(&args, "--round2");
    let r2: Value =
        serde_json::from_slice(&std::fs::read(&round2_path).unwrap()).expect("round-2 pairs");

    let rho_arms = json!([
        {"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}},
        {"name": "rho-neg-aa", "role": "control", "method": {"id": "rho.negation"}}
    ]);

    let mut written = Vec::new();
    for pair in doc["pairs"].as_array().unwrap() {
        let family = pair["family"].as_str().unwrap();
        let ell = pair["ell"].as_u64().unwrap();
        let (label, auto_arm, seed, targets) = if family == "d4" {
            (
                format!("d4-ell{}", ell),
                json!([{"name": "rho-quartic", "role": "candidate", "method": {"id": "rho.quartic"}}]),
                2026101204u64,
                32u64,
            )
        } else {
            // d19: gap sweep, rho only
            (format!("d19-ell{}", ell), Value::Null, 2026101207u64, 32u64)
        };
        let rho = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 3, rho: {} ell={}", family, ell),
            "description": "Paired planted_shared targets: the same scalars on both curves (equal r). Counts are primary; the paired quotient removes between-target variance.",
            "workloads": {
                "curves": [pair["crater"].clone(), pair["floor"].clone()],
                "targets_per_curve": targets, "target_seed": seed,
                "target_kind": "planted_shared"
            },
            "arms": rho_arms.clone(),
            "measurement": {"rounds": 1, "warmup": 1, "order": "alternate",
                             "seed": 20261012, "isolation_required": "L2",
                             "timeout_seconds": 900}
        });
        let path = format!("{}/{}-rho.json", out_dir, label);
        std::fs::write(&path, serde_json::to_vec_pretty(&rho).unwrap()).unwrap();
        written.push(path);
        if auto_arm.is_null() {
            continue;
        }
        let auto_spec = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 3, automorphism crater: {} ell={}", family, ell),
            "description": "The same crater as the rho spec (same target seed and shared law, so identical scalars), walked by the automorphism-class rho.",
            "workloads": {
                "curves": [pair["crater"].clone()],
                "targets_per_curve": targets, "target_seed": seed,
                "target_kind": "planted_shared"
            },
            "arms": auto_arm,
            "measurement": {"rounds": 1, "warmup": 1, "order": "alternate",
                             "seed": 20261012, "isolation_required": "L2",
                             "timeout_seconds": 900}
        });
        let path = format!("{}/{}-auto.json", out_dir, label);
        std::fs::write(&path, serde_json::to_vec_pretty(&auto_spec).unwrap()).unwrap();
        written.push(path);
    }
    for pair in r2["pairs"].as_array().unwrap() {
        let ell = pair["ell"].as_u64().unwrap();
        let label = format!("d3-ell{}", ell);
        let rho = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 3, paired rho: D=-3 ell={}", ell),
            "description": "Round-2 pair re-measured: the same planted scalars on both curves, 64 targets, for the paired quotient.",
            "workloads": {
                "curves": [pair["crater"].clone(), pair["floor"].clone()],
                "targets_per_curve": 64, "target_seed": 2026101205,
                "target_kind": "planted_shared"
            },
            "arms": rho_arms.clone(),
            "measurement": {"rounds": 1, "warmup": 1, "order": "alternate",
                             "seed": 20261012, "isolation_required": "L2",
                             "timeout_seconds": 900}
        });
        let cubic = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("Drastic gap round 3, paired cubic crater: D=-3 ell={}", ell),
            "description": "The same crater and the same scalars as the paired rho spec, walked by rho.cubic.",
            "workloads": {
                "curves": [pair["crater"].clone()],
                "targets_per_curve": 64, "target_seed": 2026101205,
                "target_kind": "planted_shared"
            },
            "arms": [{"name": "rho-cubic", "role": "candidate", "method": {"id": "rho.cubic"}}],
            "measurement": {"rounds": 1, "warmup": 1, "order": "alternate",
                             "seed": 20261012, "isolation_required": "L2",
                             "timeout_seconds": 900}
        });
        for (suffix, spec) in [("rho", rho), ("cubic", cubic)] {
            let path = format!("{}/{}-{}.json", out_dir, label, suffix);
            std::fs::write(&path, serde_json::to_vec_pretty(&spec).unwrap()).unwrap();
            written.push(path);
        }
    }
    for w in &written {
        println!("wrote {}", w);
    }
}

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|s| s == name)
        .expect("required argument");
    args.get(i + 1).expect("argument value").clone()
}
