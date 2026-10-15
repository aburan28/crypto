//! Verifies per-process meter receipts and summarizes complete 60-target child cost.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;

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

fn num(value: &Value, key: &str) -> u64 {
    value[key]
        .as_u64()
        .unwrap_or_else(|| panic!("missing integer {key}"))
}

fn median(mut numbers: [f64; 3]) -> f64 {
    numbers.sort_by(|a, b| a.total_cmp(b));
    numbers[1]
}

fn graph(arms: &serde_json::Map<String, Value>) -> String {
    let mut svg = String::new();
    let colors = [
        "#15803d", "#7c3aed", "#d97706", "#dc2626", "#1d4ed8", "#0891b2", "#94a3b8",
    ];
    writeln!(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" height=\"680\" viewBox=\"0 0 1200 680\">").unwrap();
    writeln!(svg, "<rect width=\"1200\" height=\"680\" fill=\"#fff\"/>").unwrap();
    writeln!(svg, "<text x=\"50\" y=\"45\" font-family=\"sans-serif\" font-size=\"28\" font-weight=\"bold\" fill=\"#102a43\">Complete 60-target solver process CPU</text>").unwrap();
    writeln!(svg, "<text x=\"50\" y=\"75\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#475569\">Three-round medians · child user + system CPU · startup and final stdout included · shared L0</text>").unwrap();
    for tick in [0, 5000, 10000, 15000, 20000] {
        let x = 280.0 + tick as f64 / 20_000.0 * 700.0;
        writeln!(
            svg,
            "<line x1=\"{x:.1}\" x2=\"{x:.1}\" y1=\"120\" y2=\"550\" stroke=\"#e2e8f0\"/>"
        )
        .unwrap();
        writeln!(svg, "<text x=\"{x:.1}\" y=\"582\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15\" fill=\"#64748b\">{tick}</text>").unwrap();
    }
    writeln!(svg, "<text x=\"630\" y=\"620\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#334155\">Child CPU (µs)</text>").unwrap();
    for (index, arm) in ARMS.iter().enumerate() {
        let row = &arms[*arm];
        let child = row["median_child_cpu_us"].as_f64().unwrap();
        let native = row["median_native_timer_us"].as_f64().unwrap();
        let y = 150 + index * 58;
        let child_width = child / 20_000.0 * 700.0;
        let native_x = 280.0 + native / 20_000.0 * 700.0;
        writeln!(svg, "<text x=\"265\" y=\"{}\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"16\" fill=\"#334155\">{arm}</text>", y + 18).unwrap();
        writeln!(svg, "<rect x=\"280\" y=\"{y}\" width=\"{child_width:.1}\" height=\"30\" rx=\"4\" fill=\"{}\" opacity=\"0.87\"/>", colors[index]).unwrap();
        writeln!(svg, "<line x1=\"{native_x:.1}\" x2=\"{native_x:.1}\" y1=\"{}\" y2=\"{}\" stroke=\"#0f172a\" stroke-width=\"3\"/>", y - 5, y + 35).unwrap();
        writeln!(svg, "<text x=\"{:.1}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"15\" fill=\"#334155\">{child:.0}</text>", 290.0 + child_width, y + 21).unwrap();
    }
    writeln!(svg, "<text x=\"50\" y=\"655\" font-family=\"sans-serif\" font-size=\"15\" fill=\"#475569\">Bar: whole child CPU. Black mark: native in-process timer in the same V4 run. Portable counted S and isolation remain open.</text>").unwrap();
    svg.push_str("</svg>\n");
    svg
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 8,
        "usage: different_j_process_summary RUN METERS REPLAY RAW_DIR SUMMARY GRAPH.svg RECEIPT.json");
    let raw = fs::read(&args[1]).unwrap();
    let meter_bytes = fs::read(&args[2]).unwrap();
    let replay_bytes = fs::read(&args[3]).unwrap();
    let replay: Value = serde_json::from_slice(&replay_bytes).unwrap();
    assert_eq!(replay["schema"], "different_j_batch_replay/v1");
    assert_eq!(replay["records_sha256"], sha256_hex(&raw));
    assert_eq!(replay["scalar_checks"], 1260);
    assert_eq!(replay["relation_witness_checks"], 10960);
    let run_text = String::from_utf8(raw.clone()).unwrap();
    let meter_text = String::from_utf8(meter_bytes.clone()).unwrap();
    let run_lines: Vec<&str> = run_text.lines().collect();
    let meter_lines: Vec<&str> = meter_text.lines().collect();
    assert_eq!((run_lines.len(), meter_lines.len()), (21, 21));
    let mut by_key = BTreeMap::<(u64, String), (Value, Value)>::new();
    for (run_line, meter_line) in run_lines.iter().zip(&meter_lines) {
        let run: Value = serde_json::from_str(run_line).unwrap();
        let meter: Value = serde_json::from_str(meter_line).unwrap();
        let round = num(&run, "round");
        let arm = run["arm"].as_str().unwrap();
        assert_eq!(run["schema"], "different_j_batch_run/v1");
        assert_eq!(run["all_verified"], true);
        assert_eq!(run["target_count"], 60);
        assert_eq!(meter["schema"], "different_j_child_process_meter/v1");
        assert_eq!(meter["round"], round);
        assert_eq!(meter["arm"], arm);
        assert_eq!(meter["exit_code"], 0);
        assert_eq!(
            num(&meter, "child_cpu_us"),
            num(&meter, "child_user_cpu_us") + num(&meter, "child_system_cpu_us")
        );
        assert!(num(&meter, "child_peak_rss_bytes") <= 512 * 1024 * 1024);
        let stem = format!("{}/round-{round}-{arm}", args[4]);
        let stdout = fs::read(format!("{stem}.stdout")).unwrap();
        let stderr = fs::read(format!("{stem}.stderr")).unwrap();
        let individual_meter: Value =
            serde_json::from_slice(&fs::read(format!("{stem}.meter.json")).unwrap()).unwrap();
        assert_eq!(stdout, format!("{run_line}\n").as_bytes());
        assert_eq!(meter, individual_meter);
        assert_eq!(num(&meter, "stdout_bytes"), stdout.len() as u64);
        assert_eq!(num(&meter, "stderr_bytes"), stderr.len() as u64);
        assert_eq!(
            fs::read_to_string(format!("{stem}.exit")).unwrap().trim(),
            "0"
        );
        assert!(by_key
            .insert((round, arm.to_string()), (run, meter))
            .is_none());
    }
    assert_eq!(by_key.len(), 21);
    let mut arm_summaries = serde_json::Map::new();
    for arm in ARMS {
        let child = std::array::from_fn::<f64, 3, _>(|round| {
            num(&by_key[&(round as u64, arm.to_string())].1, "child_cpu_us") as f64
        });
        let wall = std::array::from_fn::<f64, 3, _>(|round| {
            num(&by_key[&(round as u64, arm.to_string())].1, "wall_ns") as f64 / 1000.0
        });
        let native = std::array::from_fn::<f64, 3, _>(|round| {
            num(&by_key[&(round as u64, arm.to_string())].0, "total_cpu_ns") as f64 / 1000.0
        });
        let rss = std::array::from_fn::<u64, 3, _>(|round| {
            num(
                &by_key[&(round as u64, arm.to_string())].1,
                "child_peak_rss_bytes",
            )
        });
        let overhead = std::array::from_fn::<f64, 3, _>(|round| child[round] - native[round]);
        arm_summaries.insert(arm.to_string(), json!({
            "child_cpu_us":child,"median_child_cpu_us":median(child),
            "wall_us":wall,"median_wall_us":median(wall),
            "native_timer_us":native,"median_native_timer_us":median(native),
            "process_over_native_cpu_us":overhead,"median_process_over_native_cpu_us":median(overhead),
            "child_peak_rss_bytes":rss,
        }));
    }
    let full = arm_summaries["full-image"]["child_cpu_us"]
        .as_array()
        .unwrap();
    let bsgs = arm_summaries["source-bsgs"]["child_cpu_us"]
        .as_array()
        .unwrap();
    let rho = arm_summaries["rho-negation"]["child_cpu_us"]
        .as_array()
        .unwrap();
    let full_bsgs =
        std::array::from_fn::<f64, 3, _>(|i| full[i].as_f64().unwrap() / bsgs[i].as_f64().unwrap());
    let full_rho =
        std::array::from_fn::<f64, 3, _>(|i| full[i].as_f64().unwrap() / rho[i].as_f64().unwrap());
    let result = json!({"schema":"different_j_process_summary/v1",
        "run_sha256":sha256_hex(&raw),"meters_sha256":sha256_hex(&meter_bytes),
        "replay_sha256":sha256_hex(&replay_bytes),"arm_rounds":21,
        "scalar_checks":1260,"relation_witness_checks":10960,
        "arms":arm_summaries,
        "full_over_bsgs":{"rounds":full_bsgs,"median":median(full_bsgs)},
        "full_over_rho":{"rounds":full_rho,"median":median(full_rho)},
        "unit":"child_user_plus_system_cpu_microseconds_for_complete_60_target_process",
        "isolation":"shared_L0","portable_counted_S":null});
    fs::write(&args[5], serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    fs::write(&args[6], graph(result["arms"].as_object().unwrap())).unwrap();
    let receipt = json!({"schema":"different_j_process_summary_receipt/v1",
        "summary_sha256":sha256_hex(&fs::read(&args[5]).unwrap()),
        "graph_sha256":sha256_hex(&fs::read(&args[6]).unwrap())});
    fs::write(&args[7], serde_json::to_vec_pretty(&receipt).unwrap()).unwrap();
    println!("{receipt}");
}
