//! Independent-host meter and arithmetic comparison for the frozen fiber panel.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

const CELLS: [(u64, u64); 3] = [(457, 421), (409, 433), (421, 457)];
const PREFIXES: [u64; 4] = [1, 4, 16, 64];
const ARMS: [&str; 7] = [
    "source",
    "isogeny-image",
    "destination-direct",
    "full-image",
    "source-bsgs",
    "rho-negation",
    "rho-negation-aa",
];

type Key = (u64, u64, String);

fn number(value: &Value, key: &str) -> u64 {
    value[key]
        .as_u64()
        .unwrap_or_else(|| panic!("missing {key}"))
}

fn key(value: &Value) -> Key {
    let prefix = number(value, "prefix_length");
    let round = number(value, "round");
    let arm = value["arm"].as_str().expect("arm").to_owned();
    assert!(PREFIXES.contains(&prefix));
    assert!(round < 3);
    assert!(ARMS.contains(&arm.as_str()));
    (prefix, round, arm)
}

fn lines(path: &Path) -> (Vec<Value>, String) {
    let bytes = fs::read(path).unwrap();
    let values = bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).unwrap())
        .collect();
    (values, sha256_hex(&bytes))
}

fn indexed(values: Vec<Value>) -> BTreeMap<Key, Value> {
    let mut map = BTreeMap::new();
    for value in values {
        assert!(map.insert(key(&value), value).is_none(), "duplicate child");
    }
    assert_eq!(map.len(), 84);
    for prefix in PREFIXES {
        for round in 0..3 {
            for arm in ARMS {
                assert!(map.contains_key(&(prefix, round, arm.to_owned())));
            }
        }
    }
    map
}

fn arithmetic(value: &Value) -> Value {
    let targets: Vec<Value> = value["target_records"]
        .as_array()
        .unwrap()
        .iter()
        .map(|target| {
            json!({
                "recovered":target["recovered"],
                "query_group_ops":target["query_group_ops"],
                "rho_group_ops":target["rho_group_ops"],
                "verify_group_ops":target["verify_group_ops"],
                "target_fp_muls":target["target_fp_muls"],
                "target_fp_invs":target["target_fp_invs"],
                "trials":target["trials"],
                "relations":target["relations"]
            })
        })
        .collect();
    json!({
        "common_group_ops":value["common_group_ops"],
        "setup_group_ops":value["setup_group_ops"],
        "common_fp_muls":value["common_fp_muls"],
        "common_fp_invs":value["common_fp_invs"],
        "setup_fp_muls":value["setup_fp_muls"],
        "setup_fp_invs":value["setup_fp_invs"],
        "total_fp_muls":value["total_fp_muls"],
        "total_fp_invs":value["total_fp_invs"],
        "relation_setup":value["relation_setup"],
        "target_records":targets
    })
}

fn median(mut values: [f64; 3]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[1]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "usage: different_j_host_performance LINUX_ARTIFACT_DIR MAC_STUDY_DIR OUT.json OUT.svg"
    );
    let linux = Path::new(&args[1]);
    let mac = Path::new(&args[2]);
    assert!(!Path::new(&args[3]).exists());
    assert!(!Path::new(&args[4]).exists());
    let host_bytes = fs::read(linux.join("host.txt")).unwrap();
    let build_bytes = fs::read(linux.join("build-sha256.txt")).unwrap();
    let mut cell_summaries = Vec::new();
    let mut all_pairs = Vec::new();
    let mut total_scalars = 0;
    let mut total_relations = 0;
    let mut checked_children = 0;
    for (p, r) in CELLS {
        let cell = linux.join(format!("{p}-{r}"));
        let (runs, run_hash) = lines(&cell.join("run.jsonl"));
        let (meters, meter_hash) = lines(&cell.join("meter.jsonl"));
        let (mac_runs, mac_hash) = lines(&mac.join(format!("BATCH_{p}_{r}_V1.run.jsonl")));
        let replay_bytes = fs::read(cell.join("replay.json")).unwrap();
        let replay: Value = serde_json::from_slice(&replay_bytes).unwrap();
        assert_eq!(replay["status"], "pass");
        assert_eq!(number(&replay, "p"), p);
        assert_eq!(number(&replay, "r"), r);
        assert_eq!(number(&replay, "arm_rounds"), 84);
        assert_eq!(number(&replay, "scalar_checks"), 1785);
        assert_eq!(replay["run_sha256"], run_hash);
        assert_eq!(replay["meters_sha256"], meter_hash);
        let run = indexed(runs);
        let meter = indexed(meters);
        let old = indexed(mac_runs);
        let mut exact_arithmetic = 0u64;
        let mut cell_peak_rss = 0u64;
        for (k, v) in &run {
            let m = &meter[k];
            assert_eq!(m["exit_code"], 0);
            assert_eq!(v["all_verified"], true);
            assert_eq!(number(v, "p"), p);
            assert_eq!(number(v, "r"), r);
            assert_eq!(v["public_sha256"], replay["public_sha256"]);
            assert_eq!(
                arithmetic(v),
                arithmetic(&old[k]),
                "arithmetic drift {p}:{r} {k:?}"
            );
            exact_arithmetic += 1;
            cell_peak_rss = cell_peak_rss.max(number(m, "child_peak_rss_bytes"));
        }
        assert_eq!(exact_arithmetic, 84);
        checked_children += exact_arithmetic;
        total_scalars += number(&replay, "scalar_checks");
        total_relations += number(&replay, "relation_witness_checks");
        let mut cell_pairs = Vec::new();
        for prefix in PREFIXES {
            let mut full_rho = [0.0; 3];
            let mut full_bsgs = [0.0; 3];
            let mut aa = [0.0; 3];
            for round in 0..3 {
                let child_cpu = |arm: &str| -> f64 {
                    number(&meter[&(prefix, round, arm.to_owned())], "child_cpu_us") as f64
                };
                let f = child_cpu("full-image");
                let rho = child_cpu("rho-negation");
                let bsgs = child_cpu("source-bsgs");
                let rho_aa = child_cpu("rho-negation-aa");
                assert!(f > 0.0 && rho > 0.0 && bsgs > 0.0 && rho_aa > 0.0);
                full_rho[round as usize] = f / rho;
                full_bsgs[round as usize] = f / bsgs;
                aa[round as usize] = rho / rho_aa;
                cell_pairs.push(json!({
                    "prefix":prefix,"round":round,
                    "full_cpu_us":f,"rho_cpu_us":rho,"bsgs_cpu_us":bsgs,
                    "rho_aa_cpu_us":rho_aa,
                    "full_rho_cpu_ratio":f/rho,
                    "full_bsgs_cpu_ratio":f/bsgs,
                    "rho_aa_cpu_ratio":rho/rho_aa
                }));
            }
            all_pairs.push(json!({
                "p":p,"r":r,"prefix":prefix,
                "median_full_rho_cpu_ratio":median(full_rho),
                "median_full_bsgs_cpu_ratio":median(full_bsgs),
                "median_rho_aa_cpu_ratio":median(aa),
                "rounds":cell_pairs.iter().filter(|v|v["prefix"]==prefix).cloned().collect::<Vec<_>>()
            }));
        }
        cell_summaries.push(json!({
            "p":p,"r":r,"run_sha256":run_hash,"meter_sha256":meter_hash,
            "replay_sha256":sha256_hex(&replay_bytes),
            "mac_run_sha256":mac_hash,"exact_arithmetic_children":exact_arithmetic,
            "scalar_checks":number(&replay,"scalar_checks"),
            "relation_witness_checks":number(&replay,"relation_witness_checks"),
            "peak_child_rss_bytes":cell_peak_rss
        }));
    }
    assert_eq!(checked_children, 252);
    assert_eq!(total_scalars, 5355);
    assert_eq!(total_relations, 37061);
    let out = json!({
        "schema":"different_j_host_performance/v1",
        "status":"independent_linux_arithmetic_replay_and_descriptive_cpu",
        "wall_admission":"refused_hosted_vm_without_physical_core_reservation",
        "host_sha256":sha256_hex(&host_bytes),
        "build_sha256_manifest":sha256_hex(&build_bytes),
        "checked_children":checked_children,"checked_scalars":total_scalars,
        "checked_relation_witnesses":total_relations,
        "cells":cell_summaries,"paired_cpu":all_pairs
    });
    fs::write(
        &args[3],
        format!("{}\n", serde_json::to_string_pretty(&out).unwrap()),
    )
    .unwrap();
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="900" height="450" viewBox="0 0 900 450"><rect width="100%" height="100%" fill="white"/><text x="55" y="38" font-size="22">Independent Linux host: descriptive N=64 CPU ratios</text><text x="55" y="65" font-size="14">Three paired rounds per cell; hosted VM without physical-core isolation</text>"##,
    );
    for (ratio, y) in [(1.0, 365.0), (10.0, 265.0), (100.0, 165.0), (1000.0, 65.0)] {
        svg.push_str(&format!(
            r##"<line x1="80" y1="{y}" x2="820" y2="{y}" stroke="#ddd"/><text x="35" y="{}" font-size="14">{ratio}</text>"##,
            y + 5.0
        ));
    }
    for (field, color) in [
        ("median_full_rho_cpu_ratio", "#d95f02"),
        ("median_full_bsgs_cpu_ratio", "#4477aa"),
        ("median_rho_aa_cpu_ratio", "#1b9e77"),
    ] {
        let mut points = String::new();
        for (i, (_, r)) in CELLS.iter().enumerate() {
            let row = all_pairs
                .iter()
                .find(|v| v["r"] == *r && v["prefix"] == 64)
                .unwrap();
            let ratio = row[field].as_f64().unwrap();
            let x = 170.0 + i as f64 * 270.0;
            let y = 365.0 - 100.0 * ratio.log10();
            points.push_str(&format!("{x},{y} "));
            svg.push_str(&format!(
                r##"<circle cx="{x}" cy="{y}" r="5" fill="{color}"/>"##
            ));
        }
        svg.push_str(&format!(
            r##"<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3"/>"##
        ));
    }
    for (i, (_, r)) in CELLS.iter().enumerate() {
        let x = 150 + i * 270;
        svg.push_str(&format!(
            r##"<text x="{x}" y="395" font-size="15">r={r}</text>"##
        ));
    }
    svg.push_str(
        r##"<text x="95" y="425" font-size="14" fill="#d95f02">full/rho</text><text x="250" y="425" font-size="14" fill="#4477aa">full/BSGS</text><text x="425" y="425" font-size="14" fill="#1b9e77">rho A/A</text></svg>"##,
    );
    fs::write(&args[4], svg).unwrap();
}
