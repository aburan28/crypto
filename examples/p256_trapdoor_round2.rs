//! Run the round-2 P-256 trapdoor probes and print the recorded
//! outcomes, hypothesis by hypothesis.
//!
//! ```bash
//! cargo run --release --example p256_trapdoor_round2            # human-readable
//! cargo run --release --example p256_trapdoor_round2 -- --json  # machine-readable
//! ```
//!
//! Exits non-zero if any probe flags an anomaly, so the example doubles
//! as a regression guard on the recorded round-2 facts.

use crypto_lib::cryptanalysis::p256_trapdoor_round2::p256_trapdoor_round2_report;

fn main() {
    let json = std::env::args().any(|a| a == "--json");
    let records = p256_trapdoor_round2_report(1 << 22);

    if json {
        let rows: Vec<serde_json::Value> = records
            .iter()
            .map(|r| {
                serde_json::json!({
                    "hypothesis": r.hypothesis,
                    "name": r.name,
                    "measured": r.measured,
                    "anomaly": r.anomaly,
                })
            })
            .collect();
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "round": 2,
                "source": "RESEARCH_P256_BACKDOOR_SPECULATION.md §4 ranks 1-3",
                "probes": rows,
            }))
            .unwrap()
        );
    } else {
        for h in [1u8, 2, 3] {
            println!("## Hypothesis {h}\n");
            for r in records.iter().filter(|r| r.hypothesis == h) {
                println!("- **{}**: {}", r.name, r.measured);
                if r.anomaly {
                    println!("  - ANOMALY FLAGGED");
                }
            }
            println!();
        }
    }

    if records.iter().any(|r| r.anomaly) {
        eprintln!("round-2 anomaly flagged");
        std::process::exit(1);
    }
}
