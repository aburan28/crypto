//! Independent cut-19 confirmation gate on four fresh seeds and inactive-case controls.

use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};

const SEEDS: [&str; 4] = ["confirm_a", "confirm_b", "confirm_c", "confirm_d"];
const PRIMARY: &str = "f5_n24_m24_d4";

fn read_json(path: &Path) -> Option<Value> {
    serde_json::from_slice(&fs::read(path).ok()?).ok()
}

fn last_json_line(path: &Path) -> Option<Value> {
    let contents = fs::read_to_string(path).ok()?;
    serde_json::from_str(contents.lines().last()?).ok()
}

fn qualified(dir: &Path) -> Option<(Value, Value)> {
    let exit = read_json(&dir.join("exit.json"))?;
    let isolation = last_json_line(&dir.join("isolation.jsonl"))?;
    let paired = read_json(&dir.join("paired.json"))?;
    let clean = exit["exit_code"] == 0
        && isolation["exit_status"] == 0
        && isolation["contended_samples"] == 0
        && isolation["left_on_reserved"]["user_threads"]
            .as_array()
            .is_some_and(Vec::is_empty)
        && paired["status"] == "complete"
        && paired["runs"]
            .as_array()
            .is_some_and(|runs| runs.len() == 22 && runs.iter().all(|run| run["status"] == "ok"));
    clean.then_some((isolation, paired))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 4 || args.len() == 5,
        "usage: f5_support_cut19_confirm_analyze ROOT THREADS REPORT [CANDIDATE_RANK_BITS]"
    );
    let root = PathBuf::from(&args[1]);
    let threads: usize = args[2].parse().unwrap();
    assert!(threads == 1 || threads == 2);
    let candidate_rank_bits: u8 = args.get(4).map_or(8, |value| value.parse().unwrap());
    assert!(candidate_rank_bits == 6 || candidate_rank_bits == 8);
    let mut all_pass = true;
    let mut rows = Vec::new();
    let mut sources: Option<Value> = None;
    for seed in SEEDS {
        let mut selected = None;
        for attempt in 1..=3 {
            let dir = root.join(seed).join(format!("attempt{attempt}"));
            if let Some((_isolation, paired)) = qualified(&dir) {
                selected = Some((attempt, paired));
                break;
            }
        }
        let Some((attempt, paired)) = selected else {
            all_pass = false;
            rows.push(json!({"seed":seed, "status":"no_qualified_attempt"}));
            continue;
        };
        let hashes = json!({
            "binary_sha256": paired["binary_sha256"],
            "benchmark_source_sha256": paired["benchmark_source_sha256"],
            "harness_source_sha256": paired["harness_source_sha256"],
            "f5_source_sha256": paired["f5_source_sha256"],
            "gf2_source_sha256": paired["gf2_source_sha256"],
        });
        if let Some(first) = &sources {
            if *first != hashes {
                all_pass = false;
            }
        } else {
            sources = Some(hashes);
        }
        let primary = &paired["summary"][PRIMARY]["wall_ms"];
        let median = primary["prior_new"]["median"].as_f64().unwrap();
        let lower = primary["prior_new"]["bootstrap_95pct"][0].as_f64().unwrap();
        let aa_min = primary["aa_min"].as_f64().unwrap();
        let reference_ops = paired["runs"][0]["cases"][PRIMARY]["reduce_word_ops"]
            .as_u64()
            .unwrap();
        let candidate_ops = paired["runs"][1]["cases"][PRIMARY]["reduce_word_ops"]
            .as_u64()
            .unwrap();
        let candidate_cut = paired["runs"][1]["cases"][PRIMARY]["support_split_inner_vars"]
            .as_u64()
            .unwrap_or(0);
        let work_ratio = candidate_ops as f64 / reference_ops as f64;
        let mut errors: Vec<String> = Vec::new();
        if paired["candidate_rank_bits_cap"].as_u64() != Some(candidate_rank_bits as u64)
            || paired["runs"][1]["rank_table_bits_cap"].as_u64() != Some(candidate_rank_bits as u64)
        {
            errors.push("candidate rank-table cap mismatch".to_owned());
        }
        if paired["host"]["os"] != "linux" || paired["host"]["arch"] != "x86_64" {
            errors.push("not a Linux x86-64 host".to_owned());
        }
        if median <= 2.0 || lower <= 2.0 {
            errors.push("primary median or bootstrap lower bound <= 2.00".to_owned());
        }
        if candidate_cut != 19 {
            errors.push("candidate did not run fixed support cut 19".to_owned());
        }
        if work_ratio > 0.35 {
            errors.push("primary counted-work gate failed".to_owned());
        }
        let mut inactive_controls = Vec::new();
        for (name, case) in paired["summary"].as_object().unwrap() {
            if name == PRIMARY {
                continue;
            }
            let wall = &case["wall_ms"];
            let m = wall["prior_new"]["median"].as_f64().unwrap();
            let noise_floor = wall["aa_min"].as_f64().unwrap();
            let inactive = paired["runs"].as_array().unwrap().iter().all(|run| {
                run["cases"][name]["support_split_attempted"] == false
                    && run["cases"][name]["rank_cert_attempted"] == false
                    && run["cases"][name]["column_cert_attempted"] == false
            });
            if !inactive {
                errors.push(format!(
                    "{name}: confirmation control took a certificate route"
                ));
            }
            inactive_controls.push(json!({"case":name, "paired_median_ratio":m,
                "aa_min":noise_floor, "route_inactive":inactive}));
        }
        if !errors.is_empty() {
            all_pass = false;
        }
        rows.push(json!({
            "seed":seed, "attempt":attempt, "status":if errors.is_empty() {"pass"} else {"fail"},
            "complete_call_median_ratio":median, "bootstrap_95pct_lower":lower,
            "aa_min":aa_min, "reference_reduce_word_ops":reference_ops,
            "candidate_reduce_word_ops":candidate_ops, "candidate_inner_variables":candidate_cut,
            "counted_work_ratio":work_ratio,
            "inactive_controls":inactive_controls,
            "errors":errors,
        }));
    }
    let report = json!({"status":if all_pass {"pass"} else {"fail"},
        "threads":threads, "candidate_rank_bits_cap":candidate_rank_bits,
        "selected_sources":sources, "rows":rows});
    fs::write(
        &args[3],
        format!("{}\n", serde_json::to_string_pretty(&report).unwrap()),
    )
    .unwrap();
    println!("{}", serde_json::to_string_pretty(&report).unwrap());
    if !all_pass {
        std::process::exit(1);
    }
}
