//! Re-derive four-order counted cost trend from independently checked runs.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

const ARMS: [&str; 7] = [
    "source",
    "isogeny-image",
    "destination-direct",
    "full-image",
    "source-bsgs",
    "rho-negation",
    "rho-negation-aa",
];
const PREFIXES: [u64; 4] = [1, 4, 16, 64];

fn u(v: &Value, key: &str) -> u64 {
    v[key].as_u64().unwrap_or_else(|| panic!("missing {key}"))
}
fn gae(v: &Value) -> u64 {
    u(v, "adds") + u(v, "doubles")
}
fn median(mut values: Vec<f64>) -> f64 {
    assert_eq!(values.len(), 3);
    values.sort_by(f64::total_cmp);
    values[1]
}

fn read_cell(run_path: &str, replay_path: &str, expected_r: u64) -> (Vec<Value>, Value) {
    let bytes = fs::read(run_path).unwrap();
    let replay_bytes = fs::read(replay_path).unwrap();
    let replay: Value = serde_json::from_slice(&replay_bytes).unwrap();
    assert_eq!(replay["status"], "pass");
    assert_eq!(u(&replay, "arm_rounds"), 84);
    assert_eq!(u(&replay, "scalar_checks"), 1785);
    assert_eq!(replay["run_sha256"], sha256_hex(&bytes));
    let mut seen = BTreeMap::new();
    let mut scalar_count = 0u64;
    let mut rows = Vec::new();
    for line in bytes.split(|&b| b == b'\n').filter(|line| !line.is_empty()) {
        let record: Value = serde_json::from_slice(line).unwrap();
        assert_eq!(record["all_verified"], true);
        let p = if expected_r == 229 {
            211
        } else {
            u(&record, "p")
        };
        let r = if expected_r == 229 {
            229
        } else {
            u(&record, "r")
        };
        assert_eq!(r, expected_r);
        assert_eq!(
            record["schema"],
            if r == 229 {
                "different_j_order_batch_run/v2"
            } else {
                "different_j_trend_batch_run/v1"
            }
        );
        let prefix = u(&record, "prefix_length");
        let round = u(&record, "round");
        let arm = record["arm"].as_str().unwrap();
        assert!(PREFIXES.contains(&prefix) && round < 3 && ARMS.contains(&arm));
        assert!(seen.insert((prefix, round, arm.to_owned()), true).is_none());
        let targets = record["target_records"].as_array().unwrap();
        assert_eq!(targets.len() as u64, prefix);
        let mut group = gae(&record["common_group_ops"]) + gae(&record["setup_group_ops"]);
        let mut relations = 0u64;
        let mut trials = 0u64;
        for target in targets {
            assert_eq!(target["verified"], true);
            assert!(target["recovered"].as_u64().is_some());
            group += gae(if arm.starts_with("rho-") {
                &target["rho_group_ops"]
            } else {
                &target["query_group_ops"]
            });
            group += gae(&target["verify_group_ops"]);
            if [
                "source",
                "isogeny-image",
                "destination-direct",
                "full-image",
            ]
            .contains(&arm)
            {
                relations += u(target, "relations");
                trials += u(target, "trials");
            }
            scalar_count += 1;
        }
        let denom = prefix as f64 * (r as f64).sqrt();
        rows.push(json!({
            "p":p,"r":r,"prefix":prefix,"round":round,"arm":arm,
            "gae":group,"S_group_per_target":group as f64/denom,
            "fp_muls":u(&record,"total_fp_muls"),
            "fp_invs":u(&record,"total_fp_invs"),
            "M_per_sqrt_r_target":u(&record,"total_fp_muls") as f64/denom,
            "I_per_sqrt_r_target":u(&record,"total_fp_invs") as f64/denom,
            "cpu_ns_l0":u(&record,"total_cpu_ns"),
            "peak_rss_bytes":u(&record,"peak_rss_bytes"),
            "relations":relations,"trials":trials,
            "support_keys":record["relation_setup"]["pair_keys"],
        }));
    }
    assert_eq!(seen.len(), 84);
    assert_eq!(scalar_count, 1785);
    let cell_p = if expected_r == 229 {
        211
    } else {
        u(&rows[0], "p")
    };
    (
        rows,
        json!({"p":cell_p,
        "r":expected_r,"run_sha256":sha256_hex(&bytes),"replay_sha256":sha256_hex(&replay_bytes),
        "scalar_checks":scalar_count,"relation_witness_checks":u(&replay,"relation_witness_checks")}),
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(),11,"usage: summary RUN229 REPLAY229 RUN421 REPLAY421 RUN433 REPLAY433 RUN457 REPLAY457 OUT.json OUT.svg");
    assert!(!Path::new(&args[9]).exists() && !Path::new(&args[10]).exists());
    let mut all = Vec::new();
    let mut inputs = Vec::new();
    for (r, offset) in [(229, 1), (421, 3), (433, 5), (457, 7)] {
        let (rows, input) = read_cell(&args[offset], &args[offset + 1], r);
        all.extend(rows);
        inputs.push(input);
    }
    assert_eq!(all.len(), 336);
    let mut index = BTreeMap::<(u64, u64, u64, String), &Value>::new();
    for row in &all {
        assert!(index
            .insert(
                (
                    u(row, "r"),
                    u(row, "prefix"),
                    u(row, "round"),
                    row["arm"].as_str().unwrap().to_owned()
                ),
                row
            )
            .is_none());
    }
    let mut summary = Vec::new();
    let mut comparisons = Vec::new();
    for r in [229, 421, 433, 457] {
        for prefix in PREFIXES {
            for arm in ARMS {
                let runs: Vec<&Value> = (0..3)
                    .map(|round| index[&(r, prefix, round, arm.to_owned())])
                    .collect();
                summary.push(json!({"r":r,"prefix":prefix,"arm":arm,
                    "median_S_group_per_target":median(runs.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).collect()),
                    "min_S_group_per_target":runs.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).reduce(f64::min).unwrap(),
                    "max_S_group_per_target":runs.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).reduce(f64::max).unwrap(),
                    "median_M_per_sqrt_r_target":median(runs.iter().map(|v|v["M_per_sqrt_r_target"].as_f64().unwrap()).collect()),
                    "median_I_per_sqrt_r_target":median(runs.iter().map(|v|v["I_per_sqrt_r_target"].as_f64().unwrap()).collect()),
                    "median_cpu_ns_l0":median(runs.iter().map(|v|u(v,"cpu_ns_l0") as f64).collect()),
                    "relations_total":runs.iter().map(|v|u(v,"relations")).sum::<u64>(),
                    "trials_total":runs.iter().map(|v|u(v,"trials")).sum::<u64>(),
                    "support_keys":runs[0]["support_keys"]}));
            }
            for control in ["source-bsgs", "rho-negation"] {
                let ratios: Vec<Value> = (0..3)
                    .map(|round| {
                        let a = index[&(r, prefix, round, "full-image".to_owned())];
                        let b = index[&(r, prefix, round, control.to_owned())];
                        json!({"round":round,"group":u(a,"gae") as f64/u(b,"gae") as f64,
                        "muls":u(a,"fp_muls") as f64/u(b,"fp_muls") as f64,
                        "invs":u(a,"fp_invs") as f64/u(b,"fp_invs") as f64})
                    })
                    .collect();
                comparisons.push(json!({"r":r,"prefix":prefix,"control":control,
                    "full_exceeds_on_each_counted_axis_each_round":ratios.iter().all(|v|["group","muls","invs"].iter().all(|&k|v[k].as_f64().unwrap()>1.0)),
                    "median_group_ratio":median(ratios.iter().map(|v|v["group"].as_f64().unwrap()).collect()),
                    "ratios":ratios}));
            }
        }
    }
    let out = json!({"schema":"different_j_trend_summary/v1",
        "status":"independently_replayed_counted_axes_with_unpriced_work",
        "unit":"ecbench.gae per target per sqrt(r); Fp multiplication and inversion axes overlap group work",
        "complete_portable_S":"unpriced_nonarithmetic_work",
        "inputs":inputs,"rows":all,"summary":summary,"paired_full_vs_control":comparisons,
        "checked_children":336,"checked_scalars":7140});
    fs::write(
        &args[9],
        format!("{}\n", serde_json::to_string_pretty(&out).unwrap()),
    )
    .unwrap();

    let mut svg=String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"950\" height=\"500\" viewBox=\"0 0 950 500\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><text x=\"55\" y=\"35\" font-size=\"23\">Four checked orders: N=64 cold counted work per target</text><text x=\"55\" y=\"60\" font-size=\"14\">S_group=(adds+doubles)/(64 sqrt(r)); medians of 3 rounds; ordinal size cells and log10 y-axis</text>");
    for (value, y) in [(1, 390), (10, 290), (100, 190), (1000, 90)] {
        svg.push_str(&format!("<line x1=\"90\" y1=\"{y}\" x2=\"870\" y2=\"{y}\" stroke=\"#ddd\"/><text x=\"45\" y=\"{}\" font-size=\"14\">{value}</text>",y+5));
    }
    let colors = [
        ("full-image", "#d95f02"),
        ("rho-negation", "#1b9e77"),
        ("source-bsgs", "#4477aa"),
    ];
    for (arm, color) in colors {
        let mut points = String::new();
        for (idx, r) in [229, 421, 433, 457].iter().enumerate() {
            let row = out["summary"]
                .as_array()
                .unwrap()
                .iter()
                .find(|v| v["r"] == *r && v["prefix"] == 64 && v["arm"] == arm)
                .unwrap();
            let value = row["median_S_group_per_target"].as_f64().unwrap();
            let x = 150 + idx * 220;
            let y = 390.0 - value.log10() * 100.0;
            if !points.is_empty() {
                points.push(' ');
            }
            points.push_str(&format!("{x},{y:.1}"));
            svg.push_str(&format!("<circle cx=\"{x}\" cy=\"{y:.1}\" r=\"5\" fill=\"{color}\"/><text x=\"{}\" y=\"{:.1}\" font-size=\"13\" fill=\"{color}\">{value:.2}</text>",x+7,y-7.0));
        }
        svg.push_str(&format!("<polyline points=\"{points}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"3\"/><text x=\"680\" y=\"{}\" font-size=\"14\" fill=\"{color}\">{arm}</text>",if arm=="full-image" {35} else if arm=="rho-negation" {55} else {75}));
    }
    for (idx, r) in [229, 421, 433, 457].iter().enumerate() {
        svg.push_str(&format!(
            "<text x=\"{}\" y=\"425\" text-anchor=\"middle\" font-size=\"16\">r={r}</text>",
            150 + idx * 220
        ));
    }
    svg.push_str("<text x=\"55\" y=\"465\" font-size=\"13\">Source: 336 checked cold children, 7,140 scalar replays. Nonarithmetic work remains unpriced; L0 host timing is separate.</text></svg>\n");
    fs::write(&args[10], svg).unwrap();
}
