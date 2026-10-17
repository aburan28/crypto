//! Extend the missing rho-per-step price threshold to the selected order-1201 cell.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

const CELLS: [(u64, u64, &str); 5] = [
    (
        211,
        229,
        "d0e4616236aeb8978f5ac1549dde373be294bdf52503c954a6f964904f68faf8",
    ),
    (
        457,
        421,
        "b1054cd665cae39d07a4f64adbb1132f4d5756b59e8b3e44e1b46ce4e32009fb",
    ),
    (
        409,
        433,
        "b167728e7c0809803f387c6db8ddad4520a5fce603774a74c7d981947f1d4632",
    ),
    (
        421,
        457,
        "d75247c817c17c0bb0966f06792fb4d58f74d24434d0e83dcba45396ec98caa0",
    ),
    (
        1237,
        1201,
        "d314764e87988b2dc12675d199612f50abc538d1782761c562c2c74be5e34d89",
    ),
];
const PREFIXES: [u64; 4] = [1, 4, 16, 64];
const ARMS: [&str; 7] = [
    "destination-direct",
    "full-image",
    "isogeny-image",
    "rho-negation",
    "rho-negation-aa",
    "source",
    "source-bsgs",
];

fn u(value: &Value, field: &str) -> u64 {
    value[field]
        .as_u64()
        .unwrap_or_else(|| panic!("missing {field}"))
}
fn g(value: &Value) -> u64 {
    u(value, "adds") + u(value, "doubles")
}
fn median(mut values: [f64; 3]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[1]
}
fn read_run(path: &Path, expect: &str, replay_path: &Path) -> BTreeMap<(u64, u64, String), Value> {
    let bytes = fs::read(path).unwrap();
    assert_eq!(sha256_hex(&bytes), expect, "run digest mismatch");
    let replay_bytes = fs::read(replay_path).unwrap();
    let replay: Value = serde_json::from_slice(&replay_bytes).unwrap();
    assert_eq!(replay["status"], "pass");
    assert_eq!(replay["run_sha256"], expect);
    assert_eq!(u(&replay, "arm_rounds"), 84);
    assert_eq!(u(&replay, "scalar_checks"), 1785);
    let mut map = BTreeMap::new();
    let mut scalars = 0u64;
    for line in bytes.split(|&c| c == b'\n').filter(|line| !line.is_empty()) {
        let row: Value = serde_json::from_slice(line).unwrap();
        assert_eq!(row["all_verified"], true);
        let prefix = u(&row, "prefix_length");
        let round = u(&row, "round");
        let arm = row["arm"].as_str().unwrap().to_owned();
        assert!(PREFIXES.contains(&prefix) && round < 3 && ARMS.contains(&arm.as_str()));
        assert_eq!(
            row["target_records"].as_array().unwrap().len() as u64,
            prefix
        );
        let mut counted = g(&row["common_group_ops"]) + g(&row["setup_group_ops"]);
        for target in row["target_records"].as_array().unwrap() {
            assert_eq!(target["verified"], true);
            assert!(target["recovered"].as_u64().is_some());
            counted += g(if arm.starts_with("rho-") {
                &target["rho_group_ops"]
            } else {
                &target["query_group_ops"]
            });
            counted += g(&target["verify_group_ops"]);
            scalars += 1;
        }
        let key = (prefix, round, arm);
        assert!(map
            .insert(key, json!({"raw":row,"group":counted}))
            .is_none());
    }
    assert_eq!(map.len(), 84);
    assert_eq!(scalars, 1785);
    map
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: different_j_price_sensitivity_five STUDY_DIR OUT.json OUT.svg"
    );
    let study = Path::new(&args[1]);
    assert!(!Path::new(&args[2]).exists() && !Path::new(&args[3]).exists());
    let trend_bytes = fs::read(study.join("order1000/TREND_FIVE_V2.json")).unwrap();
    let trend: Value = serde_json::from_slice(&trend_bytes).unwrap();
    assert_eq!(u(&trend, "checked_children"), 420);
    assert_eq!(u(&trend, "checked_scalars"), 8925);
    let mut trend_index = BTreeMap::new();
    for value in trend["rows"].as_array().unwrap() {
        let key = (
            u(value, "r"),
            u(value, "prefix"),
            u(value, "round"),
            value["arm"].as_str().unwrap().to_owned(),
        );
        assert!(trend_index.insert(key, u(value, "gae")).is_none());
    }
    assert_eq!(trend_index.len(), 420);
    let mut rows = Vec::new();
    let mut n64 = Vec::new();
    for (p, r, expected_hash) in CELLS {
        let (run_path, replay_path) = if r == 229 {
            (
                study.join("../different_j_order_batch_20261015/RUN_V1.jsonl"),
                study.join("../different_j_order_batch_20261015/REPLAY_V1.json"),
            )
        } else if r == 1201 {
            (
                study.join("order1000/BATCH_1237_1201_V2.run.jsonl"),
                study.join("order1000/BATCH_1237_1201_V2.replay.json"),
            )
        } else {
            (
                study.join(format!("BATCH_{p}_{r}_V1.run.jsonl")),
                study.join(format!("BATCH_{p}_{r}_V1.replay.json")),
            )
        };
        let records = read_run(&run_path, expected_hash, &replay_path);
        let trend_input = trend["inputs"]
            .as_array()
            .unwrap()
            .iter()
            .find(|value| u(value, "r") == r)
            .unwrap();
        assert_eq!(trend_input["run_sha256"], expected_hash);
        assert_eq!(
            trend_input["replay_sha256"],
            sha256_hex(&fs::read(&replay_path).unwrap())
        );
        for prefix in PREFIXES {
            let mut thresholds = [0.0; 3];
            for round in 0..3 {
                let full = &records[&(prefix, round, "full-image".to_owned())];
                let rho = &records[&(prefix, round, "rho-negation".to_owned())];
                let full_group = u(full, "group");
                let rho_group = u(rho, "group");
                assert_eq!(
                    full_group,
                    trend_index[&(r, prefix, round, "full-image".to_owned())]
                );
                assert_eq!(
                    rho_group,
                    trend_index[&(r, prefix, round, "rho-negation".to_owned())]
                );
                let rho_steps: u64 = rho["raw"]["target_records"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .flat_map(|target| target["attempts"].as_array().unwrap().iter())
                    .map(|attempt| u(attempt, "steps"))
                    .sum();
                assert!(rho_steps > 0);
                let group_delta = full_group as i64 - rho_group as i64;
                let threshold = group_delta as f64 / rho_steps as f64;
                thresholds[round as usize] = threshold;
                rows.push(json!({
                    "p":p,"r":r,"prefix":prefix,"round":round,
                    "full_group":full_group,"rho_group":rho_group,
                    "group_delta":group_delta,"rho_steps":rho_steps,
                    "extra_rho_gae_per_step_to_tie_if_full_omitted_zero":threshold
                }));
            }
            if prefix == 64 {
                n64.push(json!({
                    "p":p,"r":r,
                    "median_gae_per_rho_step":median(thresholds),
                    "min_gae_per_rho_step":thresholds.into_iter().reduce(f64::min).unwrap(),
                    "max_gae_per_rho_step":thresholds.into_iter().reduce(f64::max).unwrap()
                }));
            }
        }
    }
    assert_eq!(rows.len(), 60);
    assert_eq!(n64.len(), 5);
    let out = json!({
        "schema":"different_j_price_sensitivity_five/v1",
        "status":"verified_missing_price_threshold_not_complete_S",
        "assumption":"all omitted full-image work set to zero; extra omitted rho cost allocated uniformly over observed walk steps",
        "unit":"group-addition equivalents per observed rho step",
        "trend_sha256":sha256_hex(&trend_bytes),
        "rows":rows,"n64_summary":n64
    });
    fs::write(
        &args[2],
        format!("{}\n", serde_json::to_string_pretty(&out).unwrap()),
    )
    .unwrap();
    let ymax = n64
        .iter()
        .map(|v| v["max_gae_per_rho_step"].as_f64().unwrap())
        .fold(0.0, f64::max)
        * 1.2;
    assert!(ymax.is_finite() && ymax > 0.0);
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="900" height="450" viewBox="0 0 900 450"><rect width="100%" height="100%" fill="white"/><text x="55" y="38" font-size="22">Five-order missing rho price threshold</text><text x="55" y="65" font-size="14">N=64 cold children; medians and round ranges; not a measured native-operation price</text>"##,
    );
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = 365.0 - 260.0 * fraction;
        svg.push_str(&format!(
            r##"<line x1="85" y1="{y}" x2="840" y2="{y}" stroke="#ddd"/><text x="28" y="{}" font-size="13">{:.0}</text>"##,
            y + 5.0,
            ymax * fraction
        ));
    }
    for (i, value) in n64.iter().enumerate() {
        let x = 120.0 + i as f64 * 155.0;
        let y = |key: &str| 365.0 - 260.0 * value[key].as_f64().unwrap() / ymax;
        let lo = y("min_gae_per_rho_step");
        let hi = y("max_gae_per_rho_step");
        let med = y("median_gae_per_rho_step");
        svg.push_str(&format!(
            r##"<line x1="{x}" y1="{lo}" x2="{x}" y2="{hi}" stroke="#d95f02" stroke-width="3"/><circle cx="{x}" cy="{med}" r="6" fill="#d95f02"/><text x="{}" y="393" font-size="15">r={}</text>"##,
            x - 22.0,
            u(value, "r")
        ));
    }
    svg.push_str(r##"<text x="130" y="425" font-size="14">Group additions of extra rho work per observed step</text></svg>"##);
    fs::write(&args[3], svg).unwrap();
}
