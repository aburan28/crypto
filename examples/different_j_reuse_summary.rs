//! Native, replay-gated summary and SVG for the different-j 60-target reuse panel.

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
const PREFIXES: [u64; 4] = [1, 4, 16, 60];

fn number(value: &Value, key: &str) -> u64 {
    value[key]
        .as_u64()
        .unwrap_or_else(|| panic!("missing integer {key}"))
}

fn median(mut numbers: [f64; 3]) -> f64 {
    numbers.sort_by(|a, b| a.total_cmp(b));
    numbers[1]
}

fn cpu_us(row: &Value, prefix_index: usize) -> f64 {
    number(&row["prefixes"][prefix_index], "total_cpu_ns") as f64 / 1000.0
}

fn ratio(rounds: [f64; 3]) -> Value {
    json!({"rounds":rounds,"median":median(rounds)})
}

fn graph(prefixes: &[Value]) -> String {
    let mut svg = String::new();
    let colors = [
        "#15803d", "#7c3aed", "#d97706", "#dc2626", "#1d4ed8", "#0891b2", "#94a3b8",
    ];
    let xs = [155.0, 390.0, 625.0, 860.0];
    let y = |value: f64| {
        610.0 - (value.ln() - 100.0f64.ln()) / (20_000.0f64.ln() - 100.0f64.ln()) * 485.0
    };
    writeln!(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1240\" height=\"760\" viewBox=\"0 0 1240 760\">").unwrap();
    writeln!(
        svg,
        "<rect width=\"1240\" height=\"760\" fill=\"#ffffff\"/>"
    )
    .unwrap();
    writeln!(svg, "<text x=\"70\" y=\"48\" font-family=\"sans-serif\" font-size=\"28\" font-weight=\"bold\" fill=\"#102a43\">Cold complete-prefix CPU over distinct targets</text>").unwrap();
    writeln!(svg, "<text x=\"70\" y=\"80\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#475569\">Three-round medians · current-thread CPU · shared L0 · input/model validation and one setup charged</text>").unwrap();
    for tick in [
        100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10_000.0, 20_000.0,
    ] {
        let yy = y(tick);
        writeln!(
            svg,
            "<line x1=\"120\" x2=\"910\" y1=\"{yy:.1}\" y2=\"{yy:.1}\" stroke=\"#e2e8f0\"/>"
        )
        .unwrap();
        writeln!(svg, "<text x=\"106\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"14\" fill=\"#64748b\">{}</text>", yy + 5.0, tick as u64).unwrap();
    }
    writeln!(svg, "<text x=\"26\" y=\"365\" transform=\"rotate(-90 26 365)\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#334155\">CPU (µs, log scale)</text>").unwrap();
    for (index, count) in PREFIXES.iter().enumerate() {
        writeln!(svg, "<text x=\"{:.1}\" y=\"644\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#334155\">{count}</text>", xs[index]).unwrap();
    }
    writeln!(svg, "<text x=\"515\" y=\"685\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"17\" fill=\"#334155\">Distinct public targets sharing one setup</text>").unwrap();
    for (arm_index, arm) in ARMS.iter().enumerate() {
        let mut values = Vec::new();
        for entry in prefixes {
            values.push(entry["arms"][*arm]["median_cpu_us"].as_f64().unwrap());
        }
        let points = values
            .iter()
            .enumerate()
            .map(|(i, &v)| format!("{:.1},{:.1}", xs[i], y(v)))
            .collect::<Vec<_>>()
            .join(" ");
        writeln!(
            svg,
            "<polyline points=\"{points}\" fill=\"none\" stroke=\"{}\" stroke-width=\"3\"/>",
            colors[arm_index]
        )
        .unwrap();
        for (i, &v) in values.iter().enumerate() {
            writeln!(
                svg,
                "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"5\" fill=\"{}\"/>",
                xs[i],
                y(v),
                colors[arm_index]
            )
            .unwrap();
        }
        let legend_y = 185 + arm_index * 44;
        writeln!(svg, "<line x1=\"965\" x2=\"1005\" y1=\"{legend_y}\" y2=\"{legend_y}\" stroke=\"{}\" stroke-width=\"4\"/>", colors[arm_index]).unwrap();
        writeln!(svg, "<text x=\"1018\" y=\"{}\" font-family=\"sans-serif\" font-size=\"16\" fill=\"#334155\">{arm}</text>", legend_y + 6).unwrap();
    }
    writeln!(svg, "<text x=\"70\" y=\"724\" font-family=\"sans-serif\" font-size=\"15\" fill=\"#64748b\">Native timer excludes process startup and final stdout; portable counted S remains open.</text>").unwrap();
    svg.push_str("</svg>\n");
    svg
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "usage: different_j_reuse_summary RUN REPLAY SUMMARY GRAPH.svg"
    );
    let raw = fs::read(&args[1]).unwrap();
    let replay: Value = serde_json::from_slice(&fs::read(&args[2]).unwrap()).unwrap();
    assert_eq!(replay["schema"], "different_j_batch_replay/v1");
    assert_eq!(replay["records_sha256"], sha256_hex(&raw));
    assert_eq!(replay["scalar_checks"], 1260);
    assert_eq!(replay["relation_witness_checks"], 10960);
    let rows: Vec<Value> = String::from_utf8(raw.clone())
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(rows.len(), 21);
    let mut by_key = BTreeMap::new();
    for row in &rows {
        assert_eq!(row["schema"], "different_j_batch_run/v1");
        assert_eq!(row["all_verified"], true);
        assert_eq!(row["target_records"].as_array().unwrap().len(), 60);
        assert_eq!(row["prefixes"].as_array().unwrap().len(), 4);
        let round = number(row, "round");
        let arm = row["arm"].as_str().unwrap();
        assert!(round < 3 && ARMS.contains(&arm));
        assert!(by_key.insert((round, arm.to_string()), row).is_none());
    }
    assert_eq!(by_key.len(), 21);
    let mut prefix_summaries = Vec::new();
    for (prefix_index, target_count) in PREFIXES.iter().enumerate() {
        let mut arm_rows = serde_json::Map::new();
        for arm in ARMS {
            let values = std::array::from_fn(|round| {
                cpu_us(by_key[&(round as u64, arm.to_string())], prefix_index)
            });
            let median_cpu = median(values);
            arm_rows.insert(
                arm.to_string(),
                json!({"round_cpu_us":values,
                "median_cpu_us":median_cpu,
                "min_cpu_us":values.iter().copied().fold(f64::INFINITY, f64::min),
                "max_cpu_us":values.iter().copied().fold(0.0, f64::max),
                "amortized_median_cpu_us_per_target":median_cpu / *target_count as f64}),
            );
        }
        let full = std::array::from_fn::<f64, 3, _>(|round| {
            cpu_us(
                by_key[&(round as u64, "full-image".to_string())],
                prefix_index,
            )
        });
        let bsgs = std::array::from_fn::<f64, 3, _>(|round| {
            cpu_us(
                by_key[&(round as u64, "source-bsgs".to_string())],
                prefix_index,
            )
        });
        let rho = std::array::from_fn::<f64, 3, _>(|round| {
            cpu_us(
                by_key[&(round as u64, "rho-negation".to_string())],
                prefix_index,
            )
        });
        prefix_summaries.push(json!({"targets":target_count,"arms":arm_rows,
            "full_over_bsgs":ratio(std::array::from_fn(|i| full[i]/bsgs[i])),
            "full_over_rho":ratio(std::array::from_fn(|i| full[i]/rho[i])),
            "full_over_best_direct":ratio(std::array::from_fn(|i| full[i]/bsgs[i].min(rho[i])))}));
    }
    let mut phase_rows = serde_json::Map::new();
    for arm in ARMS {
        let mut common = [0.0; 3];
        let mut setup = [0.0; 3];
        let mut target = [0.0; 3];
        let mut residual = [0.0; 3];
        let mut fps = [0u64; 3];
        let mut invs = [0u64; 3];
        let mut rss = [0u64; 3];
        let mut trials = [0u64; 3];
        let mut relations = [0u64; 3];
        for round in 0..3 {
            let row = by_key[&(round as u64, arm.to_string())];
            common[round] = number(row, "common_cpu_ns") as f64 / 1000.0;
            setup[round] = number(row, "setup_cpu_ns") as f64 / 1000.0;
            let target_ns: u64 = row["target_records"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| number(r, "target_cpu_ns"))
                .sum();
            target[round] = target_ns as f64 / 1000.0;
            let accounted = number(row, "common_cpu_ns") + number(row, "setup_cpu_ns") + target_ns;
            let total = number(row, "total_cpu_ns");
            assert!(total >= accounted);
            residual[round] = (total - accounted) as f64 / 1000.0;
            fps[round] = number(row, "total_fp_muls");
            invs[round] = number(row, "total_fp_invs");
            rss[round] = number(row, "peak_rss_bytes");
            trials[round] = row["target_records"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| r["trials"].as_u64().unwrap_or(0))
                .sum();
            relations[round] = row["target_records"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| r["relations"].as_u64().unwrap_or(0))
                .sum();
        }
        phase_rows.insert(
            arm.to_string(),
            json!({"common_cpu_us":common,
            "setup_cpu_us":setup,"target_cpu_us":target,"residual_cpu_us":residual,
            "median_common_cpu_us":median(common),"median_setup_cpu_us":median(setup),
            "median_target_cpu_us":median(target),"median_residual_cpu_us":median(residual),
            "field_muls":fps,"field_invs":invs,"peak_rss_bytes":rss,
            "trials":trials,"relations":relations}),
        );
    }
    let summary = json!({"schema":"different_j_batch_summary/v1",
        "run_sha256":sha256_hex(&raw),
        "replay_sha256":sha256_hex(&fs::read(&args[2]).unwrap()),
        "replayed_scalars":1260,"replayed_relation_witnesses":10960,
        "arm_rounds":21,"distinct_targets":60,
        "prefixes":prefix_summaries,"phases":phase_rows,
        "unit":"microseconds_current_thread_cpu_including_input_validation_and_setup",
        "isolation":"shared_L0","portable_counted_S":null});
    fs::write(&args[3], serde_json::to_vec_pretty(&summary).unwrap()).unwrap();
    fs::write(&args[4], graph(summary["prefixes"].as_array().unwrap())).unwrap();
    println!(
        "{}",
        json!({"schema":"different_j_batch_summary_receipt/v1",
        "summary_sha256":sha256_hex(&fs::read(&args[3]).unwrap()),
        "graph_sha256":sha256_hex(&fs::read(&args[4]).unwrap())})
    );
}
