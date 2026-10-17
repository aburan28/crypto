//! Checked native summary and SVG for the 84 cold order-229 batch children.

use std::collections::{BTreeMap, BTreeSet};
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

fn lines(bytes: &[u8]) -> Vec<Value> {
    bytes
        .split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).unwrap())
        .collect()
}

fn median(values: &[u64]) -> u64 {
    assert_eq!(values.len(), 3);
    let mut ordered = values.to_vec();
    ordered.sort_unstable();
    ordered[1]
}

fn support(record: &Value) -> BTreeSet<u64> {
    record["relation_setup"]["support_keys"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap())
        .collect()
}

fn y(value: f64, lo: f64, hi: f64) -> f64 {
    421.0 - (value.log10() - lo) / (hi - lo) * 280.0
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(),7,"usage: summary RUN.jsonl METERS.jsonl REPLAY.json CALIBRATION.jsonl SUMMARY.json GRAPH.svg");
    assert!(
        !Path::new(&args[5]).exists() && !Path::new(&args[6]).exists(),
        "refuse overwrite"
    );
    let run_bytes = fs::read(&args[1]).unwrap();
    let meter_bytes = fs::read(&args[2]).unwrap();
    let replay_bytes = fs::read(&args[3]).unwrap();
    let calibration_bytes = fs::read(&args[4]).unwrap();
    let replay: Value = serde_json::from_slice(&replay_bytes).unwrap();
    assert_eq!(replay["schema"], "different_j_order_batch_check/v2");
    assert_eq!(replay["status"], "pass");
    assert_eq!(replay["scalar_checks"], 1785);
    assert_eq!(replay["arm_rounds"], 84);
    assert_eq!(replay["run_sha256"], sha256_hex(&run_bytes));
    assert_eq!(replay["meters_sha256"], sha256_hex(&meter_bytes));
    let records = lines(&run_bytes);
    let meters = lines(&meter_bytes);
    let calibration = lines(&calibration_bytes);
    assert_eq!(
        (records.len(), meters.len(), calibration.len()),
        (84, 84, 3)
    );
    let mut ns_per_add = Vec::new();
    for (round, row) in calibration.iter().enumerate() {
        assert_eq!(row["schema"], "different_j_order_add_calibration/v1");
        assert_eq!(row["round"], round);
        assert_eq!(row["generic_add_calls"], 100_000);
        ns_per_add.push(row["cpu_ns"].as_u64().unwrap() as f64 / 100_000.0);
    }
    ns_per_add.sort_by(f64::total_cmp);
    let calibration_median = ns_per_add[1];
    assert!(calibration_median > 0.0);
    let mut run_map = BTreeMap::<(u64, u64, String), Value>::new();
    let mut meter_map = BTreeMap::<(u64, u64, String), Value>::new();
    for record in records {
        assert_eq!(record["schema"], "different_j_order_batch_run/v2");
        assert_eq!(record["all_verified"], true);
        let key = (
            record["prefix_length"].as_u64().unwrap(),
            record["round"].as_u64().unwrap(),
            record["arm"].as_str().unwrap().to_string(),
        );
        assert!(run_map.insert(key, record).is_none());
    }
    for meter in meters {
        assert_eq!(meter["schema"], "different_j_order_prefix_meter/v1");
        assert_eq!(meter["exit_code"], 0);
        let key = (
            meter["prefix_length"].as_u64().unwrap(),
            meter["round"].as_u64().unwrap(),
            meter["arm"].as_str().unwrap().to_string(),
        );
        assert!(meter_map.insert(key, meter).is_none());
    }
    assert_eq!(run_map.len(), 84);
    assert_eq!(
        run_map.keys().collect::<Vec<_>>(),
        meter_map.keys().collect::<Vec<_>>()
    );
    let mut rows = Vec::new();
    let mut ratios = Vec::new();
    let mut chart = BTreeMap::<(u64, String), (u64, u64, u64)>::new();
    for prefix in [1u64, 4, 16, 64] {
        for arm in ARMS {
            let mut cpu = Vec::new();
            let mut rss = Vec::new();
            let mut field_muls = Vec::new();
            let mut field_invs = Vec::new();
            let mut relation_total = 0u64;
            let mut trial_total = 0u64;
            let mut pair_keys = None;
            for round in 0..3u64 {
                let key = (prefix, round, arm.to_string());
                let run = &run_map[&key];
                let meter = &meter_map[&key];
                assert_eq!(run["target_count"], prefix);
                cpu.push(meter["child_cpu_us"].as_u64().unwrap());
                rss.push(meter["child_peak_rss_bytes"].as_u64().unwrap());
                field_muls.push(run["total_fp_muls"].as_u64().unwrap());
                field_invs.push(run["total_fp_invs"].as_u64().unwrap());
                if [
                    "source",
                    "isogeny-image",
                    "destination-direct",
                    "full-image",
                ]
                .contains(&arm)
                {
                    let keys = run["relation_setup"]["pair_keys"].as_u64().unwrap();
                    if let Some(previous) = pair_keys {
                        assert_eq!(previous, keys);
                    } else {
                        pair_keys = Some(keys);
                    }
                    for target in run["target_records"].as_array().unwrap() {
                        relation_total += target["relations"].as_u64().unwrap();
                        trial_total += target["trials"].as_u64().unwrap();
                    }
                }
            }
            cpu.sort_unstable();
            let cpu_median = median(&cpu);
            chart.insert((prefix, arm.to_string()), (cpu[0], cpu_median, cpu[2]));
            rows.push(json!({"prefix":prefix,"arm":arm,"rounds":3,
                "child_cpu_us":cpu,"child_cpu_median_us":cpu_median,
                "child_cpu_min_us":cpu[0],"child_cpu_max_us":cpu[2],
                "local_S_cpu_median":cpu_median as f64*1000.0/(calibration_median*(229f64).sqrt()),
                "max_rss_bytes":rss.into_iter().max().unwrap(),
                "fp_muls_median":median(&field_muls),"fp_invs_median":median(&field_invs),
                "pair_keys":pair_keys,"relations_3round_total":relation_total,
                "trials_3round_total":trial_total}));
        }
        let paired = |other: &str| -> Vec<f64> {
            (0..3u64)
                .map(|round| {
                    let full = meter_map[&(prefix, round, "full-image".to_string())]["child_cpu_us"]
                        .as_u64()
                        .unwrap() as f64;
                    let base = meter_map[&(prefix, round, other.to_string())]["child_cpu_us"]
                        .as_u64()
                        .unwrap() as f64;
                    full / base
                })
                .collect()
        };
        let mut bsgs = paired("source-bsgs");
        let mut rho = paired("rho-negation");
        bsgs.sort_by(f64::total_cmp);
        rho.sort_by(f64::total_cmp);
        ratios.push(
            json!({"prefix":prefix,"full_over_bsgs_paired":bsgs,"full_over_rho_paired":rho,
            "full_over_bsgs_median":bsgs[1],"full_over_rho_median":rho[1]}),
        );
    }
    let direct = support(&run_map[&(64, 0, "destination-direct".to_string())]);
    let full = support(&run_map[&(64, 0, "full-image".to_string())]);
    let intersection = direct.intersection(&full).count();
    let summary = json!({"schema":"different_j_order_batch_summary/v1",
        "run_sha256":sha256_hex(&run_bytes),"meter_sha256":sha256_hex(&meter_bytes),
        "replay_sha256":sha256_hex(&replay_bytes),"calibration_sha256":sha256_hex(&calibration_bytes),
        "calibration_ns_per_add_sorted":ns_per_add,"calibration_ns_per_add_median":calibration_median,
        "rows":rows,"paired_ratios":ratios,
        "destination_support":{"direct":direct.len(),"full_image":full.len(),
            "intersection":intersection,"direct_only":direct.len()-intersection,"full_only":full.len()-intersection},
        "isolation":"shared L0","portable_counted_S":Value::Null,
        "local_S_cpu_status":"derived from measured CPU/add conversion; not a portable bound"});
    fs::write(
        &args[5],
        format!("{}\n", serde_json::to_string_pretty(&summary).unwrap()),
    )
    .unwrap();

    let plot_arms = [
        ("full-image", "#ea580c"),
        ("source-bsgs", "#0e7490"),
        ("rho-negation", "#64748b"),
    ];
    let values: Vec<f64> = chart
        .values()
        .flat_map(|&(a, b, c)| [a, b, c].map(|x| x as f64))
        .collect();
    let lo = values
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
        .log10()
        .floor();
    let hi = values.iter().copied().fold(0f64, f64::max).log10().ceil();
    assert!(hi > lo);
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="535" viewBox="0 0 1000 535"><rect width="1000" height="535" fill="#f8fafc"/><text x="54" y="36" font-family="Arial,sans-serif" font-size="22" font-weight="700" fill="#0f172a">Order-229 cold child CPU by distinct-target prefix</text><text x="54" y="59" font-family="Arial,sans-serif" font-size="13" fill="#475569">Medians and min–max of 3 arm-rotated rounds · same 64-target public order · shared L0 host</text>"##,
    );
    for exp in (lo as i32)..=(hi as i32) {
        let val = 10f64.powi(exp);
        let py = y(val, lo, hi);
        svg.push_str(&format!("<line x1=\"95\" y1=\"{py:.1}\" x2=\"955\" y2=\"{py:.1}\" stroke=\"#cbd5e1\" stroke-dasharray=\"3 4\"/><text x=\"47\" y=\"{:.1}\" font-family=\"Arial,sans-serif\" font-size=\"12\" fill=\"#475569\">10^{exp}</text>",py+4.0));
    }
    for (pindex, prefix) in [1u64, 4, 16, 64].into_iter().enumerate() {
        let cx = 200.0 + pindex as f64 * 220.0;
        svg.push_str(&format!("<text x=\"{cx:.1}\" y=\"455\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">N={prefix}</text>"));
        for (aindex, (arm, color)) in plot_arms.iter().enumerate() {
            let (min, med, max) = chart[&(prefix, arm.to_string())];
            let x = cx - 52.0 + aindex as f64 * 52.0;
            let top = y(med as f64, lo, hi);
            let bottom = y(10f64.powf(lo), lo, hi);
            let upper = y(max as f64, lo, hi);
            let lower = y(min as f64, lo, hi);
            svg.push_str(&format!("<rect x=\"{x:.1}\" y=\"{top:.1}\" width=\"30\" height=\"{:.1}\" fill=\"{color}\" opacity=\"0.9\" rx=\"2\"/><line x1=\"{:.1}\" y1=\"{upper:.1}\" x2=\"{:.1}\" y2=\"{lower:.1}\" stroke=\"#0f172a\" stroke-width=\"1.5\"/><line x1=\"{:.1}\" y1=\"{upper:.1}\" x2=\"{:.1}\" y2=\"{upper:.1}\" stroke=\"#0f172a\"/><line x1=\"{:.1}\" y1=\"{lower:.1}\" x2=\"{:.1}\" y2=\"{lower:.1}\" stroke=\"#0f172a\"/><text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"10\" fill=\"#334155\">{med}</text>",bottom-top,x+15.0,x+15.0,x+9.0,x+21.0,x+9.0,x+21.0,x+15.0,upper-7.0));
        }
    }
    svg.push_str(r##"<text x="30" y="112" font-family="Arial,sans-serif" font-size="12" fill="#475569">child CPU (µs, log scale)</text><rect x="110" y="484" width="14" height="14" fill="#ea580c"/><text x="130" y="496" font-family="Arial,sans-serif" font-size="13" fill="#334155">full rational image</text><rect x="328" y="484" width="14" height="14" fill="#0e7490"/><text x="348" y="496" font-family="Arial,sans-serif" font-size="13" fill="#334155">reusable source BSGS</text><rect x="565" y="484" width="14" height="14" fill="#64748b"/><text x="585" y="496" font-family="Arial,sans-serif" font-size="13" fill="#334155">source negation rho</text><text x="110" y="522" font-family="Arial,sans-serif" font-size="12" fill="#64748b">Whole-child CPU includes construction, setup, queries, recovery, verification, startup and output; L0 descriptive.</text></svg>"##);
    fs::write(&args[6], svg).unwrap();
    println!(
        "{}",
        json!({"status":"ok","arm_rounds":84,"scalar_checks":1785,
        "support_intersection":intersection,"calibration_ns_per_add_median":calibration_median})
    );
}
