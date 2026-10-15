//! Re-derive the portable counted axes from the independently checked V1 run.
//! Protocol: research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/PROTOCOL_V1.md.

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
const R: u64 = 229;

fn u(v: &Value, name: &str) -> u64 {
    v[name].as_u64().unwrap_or_else(|| panic!("missing {name}"))
}

fn gae(v: &Value) -> u64 {
    u(v, "adds") + u(v, "doubles")
}

fn median(mut values: Vec<f64>) -> f64 {
    assert_eq!(values.len(), 3);
    values.sort_by(f64::total_cmp);
    values[1]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "usage: different_j_order_portable RUN.jsonl REPLAY.json OUT.json OUT.svg"
    );
    assert!(!Path::new(&args[3]).exists() && !Path::new(&args[4]).exists());
    let run_bytes = fs::read(&args[1]).unwrap();
    let replay: Value = serde_json::from_slice(&fs::read(&args[2]).unwrap()).unwrap();
    assert_eq!(replay["schema"], "different_j_order_batch_check/v2");
    assert_eq!(replay["status"], "pass");
    assert_eq!(u(&replay, "arm_rounds"), 84);
    assert_eq!(u(&replay, "scalar_checks"), 1785);
    assert_eq!(replay["run_sha256"], sha256_hex(&run_bytes));

    let mut rows = BTreeMap::<(u64, u64, String), Value>::new();
    let root = (R as f64).sqrt();
    let mut scalar_count = 0u64;
    for line in run_bytes
        .split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
    {
        let record: Value = serde_json::from_slice(line).unwrap();
        assert_eq!(record["schema"], "different_j_order_batch_run/v2");
        assert_eq!(record["all_verified"], true);
        let prefix = u(&record, "prefix_length");
        let round = u(&record, "round");
        let arm = record["arm"].as_str().unwrap().to_owned();
        assert!(PREFIXES.contains(&prefix) && round < 3 && ARMS.contains(&arm.as_str()));
        let targets = record["target_records"].as_array().unwrap();
        assert_eq!(targets.len() as u64, prefix);
        let mut counted = gae(&record["common_group_ops"]) + gae(&record["setup_group_ops"]);
        for target in targets {
            assert_eq!(target["verified"], true);
            assert!(target["recovered"].as_u64().is_some());
            let work = if arm.starts_with("rho-") {
                &target["rho_group_ops"]
            } else {
                &target["query_group_ops"]
            };
            counted += gae(work) + gae(&target["verify_group_ops"]);
            scalar_count += 1;
        }
        let m = u(&record, "total_fp_muls");
        let i = u(&record, "total_fp_invs");
        let key = (prefix, round, arm.clone());
        let row = json!({
            "prefix":prefix,"round":round,"arm":arm,"r":R,
            "targets":prefix,"gae":counted,"fp_muls":m,"fp_invs":i,
            "S_group_per_target":counted as f64/(prefix as f64*root),
            "M_per_sqrt_r_target":m as f64/(prefix as f64*root),
            "I_per_sqrt_r_target":i as f64/(prefix as f64*root),
            "peak_rss_bytes":u(&record,"peak_rss_bytes"),
            "unpriced":["hashes","table_lookup_insert","row_arithmetic","allocation","startup","input_output","field_additions"],
        });
        assert!(rows.insert(key, row).is_none(), "duplicate child");
    }
    assert_eq!(rows.len(), 84);
    assert_eq!(scalar_count, 1785);
    let mut summary = Vec::new();
    let mut pair_ratios = Vec::new();
    for prefix in PREFIXES {
        for arm in ARMS {
            let row_set: Vec<&Value> = (0..3)
                .map(|round| &rows[&(prefix, round, arm.to_owned())])
                .collect();
            summary.push(json!({
                "prefix":prefix,"arm":arm,
                "median_S_group_per_target":median(row_set.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).collect()),
                "min_S_group_per_target":row_set.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).reduce(f64::min).unwrap(),
                "max_S_group_per_target":row_set.iter().map(|v|v["S_group_per_target"].as_f64().unwrap()).reduce(f64::max).unwrap(),
                "median_M_per_sqrt_r_target":median(row_set.iter().map(|v|v["M_per_sqrt_r_target"].as_f64().unwrap()).collect()),
                "median_I_per_sqrt_r_target":median(row_set.iter().map(|v|v["I_per_sqrt_r_target"].as_f64().unwrap()).collect()),
            }));
        }
        for control in ["source-bsgs", "rho-negation"] {
            let ratios: Vec<Value> = (0..3)
                .map(|round| {
                    let full = &rows[&(prefix, round, "full-image".to_owned())];
                    let base = &rows[&(prefix, round, control.to_owned())];
                    json!({"round":round,
                        "group":u(full,"gae") as f64/u(base,"gae") as f64,
                        "muls":u(full,"fp_muls") as f64/u(base,"fp_muls") as f64,
                        "invs":u(full,"fp_invs") as f64/u(base,"fp_invs") as f64})
                })
                .collect();
            let dominates = ratios.iter().all(|v| {
                ["group", "muls", "invs"]
                    .iter()
                    .all(|&k| v[k].as_f64().unwrap() > 1.0)
            });
            pair_ratios.push(json!({"prefix":prefix,"control":control,"full_exceeds_control_each_axis_each_round":dominates,"ratios":ratios}));
        }
    }
    let out = json!({
        "schema":"different_j_order_portable/v1",
        "status":"checked_counted_axes_with_unpriced_work",
        "r":R,"floor_per_target":std::f64::consts::PI.sqrt()/2.0,
        "run_sha256":sha256_hex(&run_bytes),"replay_sha256":sha256_hex(&fs::read(&args[2]).unwrap()),
        "checked_children":rows.len(),"checked_scalars":scalar_count,
        "unit":"ecbench.gae per target per sqrt(r); separate Fp multiplication and inversion axes overlap group work",
        "rows":rows.values().collect::<Vec<_>>(),"summary":summary,"paired_full_vs_control":pair_ratios,
        "complete_portable_S":"unpriced_nonarithmetic_work",
    });
    fs::write(
        &args[3],
        format!("{}\n", serde_json::to_string_pretty(&out).unwrap()),
    )
    .unwrap();

    let mut svg = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1050\" height=\"600\" viewBox=\"0 0 1050 600\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><text x=\"45\" y=\"35\" font-size=\"23\">Order 229: cold counted group work per target</text><text x=\"45\" y=\"60\" font-size=\"14\">S = additions + doublings divided by N sqrt(229); median of 3 rounds; log scale</text>");
    let n64: Vec<&Value> = out["summary"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|v| v["prefix"] == 64)
        .collect();
    for (idx, row) in n64.iter().enumerate() {
        let value = row["median_S_group_per_target"].as_f64().unwrap();
        let width = ((value.log10() + 1.0) / 5.0 * 620.0).clamp(0.0, 620.0);
        let y = 95 + idx * 68;
        svg.push_str(&format!("<text x=\"45\" y=\"{}\" font-size=\"16\">{}</text><rect x=\"280\" y=\"{}\" width=\"{:.1}\" height=\"30\" fill=\"{}\"/><text x=\"{}\" y=\"{}\" font-size=\"15\">{:.2}</text>", y+21,row["arm"].as_str().unwrap(),y,width,if row["arm"]=="full-image" {"#d95f02"} else {"#1b9e77"},290.0+width,y+21,value));
    }
    svg.push_str("<text x=\"45\" y=\"580\" font-size=\"12\">Source: checked V1 84-child ledger; field and nonarithmetic work remain separate. This is not complete S.</text></svg>\n");
    fs::write(&args[4], svg).unwrap();
}
