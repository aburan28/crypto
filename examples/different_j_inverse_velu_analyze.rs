//! Native paired analysis of the frozen direct inverse-Velu 60-child panel.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

const PREFIXES: [u64; 4] = [1, 4, 16, 64];
const ARMS: [&str; 5] = [
    "inverse-velu",
    "full-image",
    "source-bsgs",
    "rho-negation",
    "rho-negation-aa",
];
type Key = (u64, u64, String);

fn number(v: &Value, key: &str) -> u64 {
    v[key].as_u64().unwrap_or_else(|| panic!("missing {key}"))
}

fn key(v: &Value) -> Key {
    let prefix = number(v, "prefix_length");
    let round = number(v, "round");
    let arm = v["arm"].as_str().unwrap().to_owned();
    assert!(PREFIXES.contains(&prefix) && round < 3 && ARMS.contains(&arm.as_str()));
    (prefix, round, arm)
}

fn indexed(path: &Path) -> (BTreeMap<Key, Value>, String) {
    let bytes = fs::read(path).unwrap();
    let mut out = BTreeMap::new();
    for line in bytes
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
    {
        let value: Value = serde_json::from_slice(line).unwrap();
        assert!(out.insert(key(&value), value).is_none(), "duplicate child");
    }
    assert_eq!(out.len(), 60);
    for prefix in PREFIXES {
        for round in 0..3 {
            for arm in ARMS {
                assert!(out.contains_key(&(prefix, round, arm.to_owned())));
            }
        }
    }
    (out, sha256_hex(&bytes))
}

fn group_ops(v: &Value) -> u64 {
    if v.is_null() {
        0
    } else {
        number(v, "adds") + number(v, "doubles")
    }
}

fn counted_group_total(v: &Value) -> u64 {
    let mut total = group_ops(&v["common_group_ops"]) + group_ops(&v["setup_group_ops"]);
    for target in v["target_records"].as_array().unwrap() {
        total += group_ops(&target["query_group_ops"]);
        total += group_ops(&target["rho_group_ops"]);
        total += group_ops(&target["verify_group_ops"]);
    }
    total
}

fn median(mut values: Vec<f64>) -> f64 {
    assert_eq!(values.len(), 3);
    values.sort_by(f64::total_cmp);
    values[1]
}

fn range(values: &[f64]) -> [f64; 2] {
    [
        values.iter().copied().fold(f64::INFINITY, f64::min),
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    ]
}

fn metric(v: &Value, m: &Value, prefix: u64, round: u64) -> Value {
    assert_eq!(v["all_verified"], true);
    assert_eq!(m["exit_code"], 0);
    assert_eq!(v["p"], 1237);
    assert_eq!(v["r"], 1201);
    assert_eq!(v["target_count"], prefix);
    assert_eq!(v["round"], round);
    let group = counted_group_total(v);
    let cpu = number(m, "child_cpu_us");
    let fp_muls = number(v, "total_fp_muls");
    let fp_invs = number(v, "total_fp_invs");
    assert!(cpu > 0 && group > 0);
    json!({
        "prefix":prefix,"round":round,
        "group_ops":group,
        "s_group":group as f64/(prefix as f64*(1201f64).sqrt()),
        "fp_muls":fp_muls,"fp_invs":fp_invs,
        "child_cpu_us":cpu,
        "child_peak_rss_bytes":number(m,"child_peak_rss_bytes"),
    })
}

fn paired_ratio(
    metrics: &BTreeMap<Key, Value>,
    prefix: u64,
    numerator: &str,
    denominator: &str,
    field: &str,
) -> Value {
    let mut values = Vec::new();
    for round in 0..3 {
        let a = metrics[&(prefix, round, numerator.to_owned())][field]
            .as_f64()
            .unwrap();
        let b = metrics[&(prefix, round, denominator.to_owned())][field]
            .as_f64()
            .unwrap();
        assert!(a > 0.0 && b > 0.0);
        values.push(a / b);
    }
    json!({
        "numerator":numerator,"denominator":denominator,"metric":field,
        "prefix":prefix,"round_ratios":values,
        "median_ratio":median(values.clone()),"round_range":range(&values)
    })
}

fn svg(pairs: &[Value]) -> String {
    let wanted = [
        ("full-image", "inverse-velu", "full / inverse", "#d95f02"),
        ("inverse-velu", "rho-negation", "inverse / rho", "#4477aa"),
        ("inverse-velu", "source-bsgs", "inverse / BSGS", "#aa4499"),
        (
            "rho-negation",
            "rho-negation-aa",
            "rho / rho A/A",
            "#1b9e77",
        ),
    ];
    let mut out = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="510" viewBox="0 0 1050 510"><rect width="100%" height="100%" fill="white"/><text x="45" y="38" font-size="22">Inverse-Velu: cold N=64 process CPU ratios</text><text x="45" y="65" font-size="14">Three paired rounds; Darwin shared-host L0; all targets replayed</text>"##,
    );
    let y_of = |ratio: f64| 405.0 - 80.0 * ratio.log10();
    for (ratio, label) in [(1.0, "1"), (10.0, "10"), (100.0, "100")] {
        let y = y_of(ratio);
        out.push_str(&format!(r##"<line x1="85" y1="{y}" x2="1000" y2="{y}" stroke="#ddd"/><text x="40" y="{}" font-size="14">{label}</text>"##,y+5.0));
    }
    for (index, (num, den, label, color)) in wanted.iter().enumerate() {
        let row = pairs
            .iter()
            .find(|v| {
                v["prefix"] == 64
                    && v["metric"] == "child_cpu_us"
                    && v["numerator"] == *num
                    && v["denominator"] == *den
            })
            .unwrap();
        let mid = row["median_ratio"].as_f64().unwrap();
        let low = row["round_range"][0].as_f64().unwrap();
        let high = row["round_range"][1].as_f64().unwrap();
        assert!(low > 0.0 && high < 100.0);
        let x = 190.0 + 245.0 * index as f64;
        out.push_str(&format!(r##"<line x1="{x}" y1="{}" x2="{x}" y2="{}" stroke="{color}" stroke-width="3"/><line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{color}" stroke-width="3"/><line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{color}" stroke-width="3"/><circle cx="{x}" cy="{}" r="7" fill="{color}"/><text x="{}" y="445" font-size="14">{label}</text><text x="{}" y="469" font-size="13" fill="{color}">median {:.3}x</text>"##,
            y_of(low),y_of(high),x-12.0,y_of(low),x+12.0,y_of(low),x-12.0,y_of(high),x+12.0,y_of(high),y_of(mid),x-65.0,x-42.0,mid));
    }
    out.push_str("</svg>");
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "usage: different_j_inverse_velu_analyze RUN_DIR OUT.json OUT.svg"
    );
    let run_dir = Path::new(&args[1]);
    assert!(!Path::new(&args[2]).exists() && !Path::new(&args[3]).exists());
    let (runs, run_hash) = indexed(&run_dir.join("run.jsonl"));
    let (meters, meter_hash) = indexed(&run_dir.join("meter.jsonl"));
    let replay: Value =
        serde_json::from_slice(&fs::read(run_dir.join("replay.json")).unwrap()).unwrap();
    assert_eq!(replay["status"], "pass");
    assert_eq!(replay["arm_rounds"], 60);
    assert_eq!(replay["scalar_checks"], 1275);
    assert_eq!(replay["run_sha256"], run_hash);
    assert_eq!(replay["meters_sha256"], meter_hash);
    let mut metrics = BTreeMap::new();
    for (key, run) in &runs {
        assert_eq!(run["construction"], "inverse-velu-v2");
        let value = metric(run, &meters[key], key.0, key.1);
        assert!(metrics.insert(key.clone(), value).is_none());
    }
    let mut matched_relation_rows = 0u64;
    for prefix in PREFIXES {
        for round in 0..3 {
            let a = &runs[&(prefix, round, "inverse-velu".to_owned())];
            let b = &runs[&(prefix, round, "full-image".to_owned())];
            assert_eq!(
                a["relation_setup"]["support_keys"],
                b["relation_setup"]["support_keys"]
            );
            assert_eq!(
                a["relation_setup"]["pair_entries"],
                b["relation_setup"]["pair_entries"]
            );
            for (left, right) in a["target_records"]
                .as_array()
                .unwrap()
                .iter()
                .zip(b["target_records"].as_array().unwrap())
            {
                for field in [
                    "recovered",
                    "trials",
                    "relations",
                    "rank",
                    "relation_witnesses",
                ] {
                    assert_eq!(
                        left[field], right[field],
                        "relation-row drift {prefix} {round} {field}"
                    );
                }
                matched_relation_rows += 1;
            }
        }
    }
    assert_eq!(matched_relation_rows, 255);
    let mut rows = Vec::new();
    for prefix in PREFIXES {
        for arm in ARMS {
            let values: Vec<Value> = (0..3)
                .map(|round| metrics[&(prefix, round, arm.to_owned())].clone())
                .collect();
            let median_field =
                |field: &str| median(values.iter().map(|v| v[field].as_f64().unwrap()).collect());
            rows.push(json!({"prefix":prefix,"arm":arm,
                "median_s_group":median_field("s_group"),
                "median_group_ops":median_field("group_ops"),
                "median_fp_muls":median_field("fp_muls"),
                "median_fp_invs":median_field("fp_invs"),
                "median_child_cpu_us":median_field("child_cpu_us"),
                "max_child_peak_rss_bytes":values.iter().map(|v|number(v,"child_peak_rss_bytes")).max().unwrap(),
                "rounds":values}));
        }
    }
    let mut pairs = Vec::new();
    for prefix in PREFIXES {
        for field in ["group_ops", "fp_muls", "fp_invs", "child_cpu_us"] {
            for (num, den) in [
                ("full-image", "inverse-velu"),
                ("inverse-velu", "rho-negation"),
                ("inverse-velu", "source-bsgs"),
                ("rho-negation", "rho-negation-aa"),
            ] {
                pairs.push(paired_ratio(&metrics, prefix, num, den, field));
            }
        }
    }
    let result = json!({"schema":"different_j_inverse_velu_analysis/v1",
        "status":"complete_replay_with_unpriced_native_categories",
        "host_level":"L0_shared_Darwin_descriptive_cpu",
        "portable_s_status":"lower_bound_only",
        "run_sha256":run_hash,"meter_sha256":meter_hash,
        "checker_sha256":sha256_hex(&fs::read(run_dir.join("replay.json")).unwrap()),
        "checked_children":60,"checked_scalars":1275,
        "checked_relation_witnesses":number(&replay,"relation_witness_checks"),
        "matched_inverse_full_relation_rows":matched_relation_rows,
        "rows":rows,"paired_ratios":pairs});
    fs::write(
        &args[2],
        format!("{}\n", serde_json::to_string_pretty(&result).unwrap()),
    )
    .unwrap();
    fs::write(&args[3], svg(result["paired_ratios"].as_array().unwrap())).unwrap();
}
