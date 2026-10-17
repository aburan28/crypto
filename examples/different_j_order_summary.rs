//! Native summary and SVG from the replayed larger-order different-j screen.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 5, "raw checker summary graph paths");
    let raw = fs::read(&args[1]).expect("raw screen");
    let replay = fs::read(&args[2]).expect("independent replay");
    let replay_json: Value = serde_json::from_slice(&replay).unwrap();
    assert_eq!(replay_json["status"], "pass");
    assert!(!Path::new(&args[3]).exists(), "refuse summary overwrite");
    assert!(!Path::new(&args[4]).exists(), "refuse graph overwrite");
    let records: Vec<Value> = raw
        .split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| serde_json::from_slice(line).unwrap())
        .collect();
    let final_row = records.last().expect("summary row");
    assert_eq!(final_row["schema"], "different_j_order_screen_summary/v1");
    let rows = &records[..records.len() - 1];
    assert_eq!(
        rows.len(),
        replay_json["candidates"].as_u64().unwrap() as usize
    );
    let mut reason_counts = BTreeMap::<String, usize>::new();
    let mut chart_rows = Vec::new();
    let mut accepted = Vec::new();
    for row in rows {
        let p = row["p"].as_u64().unwrap();
        let status = row["status"].as_str().unwrap();
        let reason = row["reason"].as_str().unwrap_or("accept");
        *reason_counts.entry(reason.to_string()).or_default() += 1;
        let max_factor = row["factors"]
            .as_array()
            .and_then(|a| a.iter().filter_map(|x| x[0].as_u64()).max());
        if status == "accept" {
            accepted.push(json!({"p":p,"r":row["r"],"n1":row["n1"],
                "n2":row["n2"],"cofactor":row["cofactor"],
                "source_j":row["source"]["j"],"destination_j":row["destination"]["j"]}));
        }
        chart_rows.push(json!({"p":p,"status":status,"reason":reason,
            "largest_prime_factor":max_factor,"r":row["r"]}));
    }
    assert_eq!(
        accepted.len(),
        replay_json["accepted"].as_u64().unwrap() as usize
    );
    assert_eq!(final_row["selected"], replay_json["selected"]);
    let summary = json!({
        "schema":"different_j_order_summary/v1",
        "raw_sha256":sha256_hex(&raw),"replay_sha256":sha256_hex(&replay),
        "candidate_count":rows.len(),"accepted_count":accepted.len(),
        "map_witnesses":replay_json["map_witnesses"],
        "selected":final_row["selected"],"accepted":accepted,
        "reason_counts":reason_counts,"rows":chart_rows,
        "stage":"bounded_candidate_screen","portable_counted_S":Value::Null
    });
    fs::write(
        &args[3],
        format!("{}\n", serde_json::to_string_pretty(&summary).unwrap()),
    )
    .unwrap();

    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="430" viewBox="0 0 1000 430">
<rect width="1000" height="430" fill="#f8fafc"/>
<text x="54" y="34" font-family="Arial,sans-serif" font-size="22" font-weight="700" fill="#0f172a">Different-j order screen: p ≤ 271</text>
<text x="54" y="56" font-family="Arial,sans-serif" font-size="13" fill="#475569">Largest prime factor of #E(Fₚ²), 20 eligible primes; green bars pass r ≥ 211 with multiplicity one</text>
"##,
    );
    let plot_y = 352.0f64;
    let scale = 270.0 / 260.0;
    for y in [0u64, 100, 211, 260] {
        let py = plot_y - y as f64 * scale;
        let dash = if y == 211 { "6 4" } else { "2 5" };
        svg.push_str(&format!(
            "<line x1=\"54\" y1=\"{py:.1}\" x2=\"968\" y2=\"{py:.1}\" stroke=\"#94a3b8\" stroke-dasharray=\"{dash}\"/><text x=\"12\" y=\"{:.1}\" font-family=\"Arial,sans-serif\" font-size=\"12\" fill=\"#475569\">{y}</text>\n",
            py + 4.0
        ));
    }
    for (i, row) in chart_rows.iter().enumerate() {
        let p = row["p"].as_u64().unwrap();
        let factor = row["largest_prime_factor"].as_u64().unwrap_or(0);
        let x = 62.0 + i as f64 * 45.0;
        let height = factor as f64 * scale;
        let y = plot_y - height;
        let color = if p == 211 {
            "#f97316"
        } else if row["status"] == "accept" {
            "#16a34a"
        } else {
            "#64748b"
        };
        svg.push_str(&format!(
            "<rect x=\"{x:.1}\" y=\"{y:.1}\" width=\"31\" height=\"{height:.1}\" fill=\"{color}\" rx=\"2\"/><text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"11\" font-weight=\"600\" fill=\"#0f172a\">{factor}</text><text x=\"{:.1}\" y=\"374\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"11\" fill=\"#334155\">{p}</text>\n",
            x + 15.5,
            y - 5.0,
            x + 15.5
        ));
    }
    svg.push_str(
        r##"<rect x="54" y="393" width="13" height="13" fill="#f97316"/><text x="74" y="404" font-family="Arial,sans-serif" font-size="12" fill="#334155">selected p=211, r=229</text>
<rect x="270" y="393" width="13" height="13" fill="#16a34a"/><text x="290" y="404" font-family="Arial,sans-serif" font-size="12" fill="#334155">other qualifying pair</text>
<rect x="473" y="393" width="13" height="13" fill="#64748b"/><text x="493" y="404" font-family="Arial,sans-serif" font-size="12" fill="#334155">subgroup rejection</text>
<text x="707" y="404" font-family="Arial,sans-serif" font-size="12" fill="#475569">candidate census, not DLP cost</text></svg>
"##,
    );
    fs::write(&args[4], svg).unwrap();
    println!(
        "{}",
        json!({"status":"ok","candidates":rows.len(),
        "accepted":accepted.len(),"selected":final_row["selected"]})
    );
}
