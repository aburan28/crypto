//! Native reduction of the frozen seven-arm held-out JSON and meter records.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use serde_json::{json, Value};

const ARMS: [&str; 9] = [
    "chi-resultant", "direct-pair", "source-s3", "quotient-s3",
    "inverse-velu", "full-image", "source-bsgs", "rho-negation", "rho-negation-aa",
];

fn number(row: &Value, key: &str) -> u64 {
    row[key].as_u64().unwrap_or_else(|| panic!("missing {key}"))
}

fn group_count(row: &Value) -> u64 {
    ["setup_group_ops", "query_group_ops", "verify_group_ops"]
        .iter()
        .map(|phase| number(&row[phase], "adds") + number(&row[phase], "doubles"))
        .sum()
}

#[derive(Default)]
struct Totals {
    count: u64,
    group: u64,
    meter_cpu_us: u64,
    meter_wall_ns: u64,
    solver_cpu_ns: u64,
    peak_rss_bytes: u64,
    trials: u64,
    relations: u64,
    rank: u64,
    fp_muls: u64,
    fp_invs: u64,
}

fn mean(total: u64, count: u64) -> f64 {
    total as f64 / count as f64
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 2, "usage: summary RUN_DIRECTORY");
    let root = Path::new(&args[1]);
    let mut totals: BTreeMap<String, Totals> = BTreeMap::new();
    let mut paired: BTreeMap<(u64, u64), BTreeMap<String, (u64, u64, u64, u64)>> = BTreeMap::new();
    for round in 0..3u64 {
        for arm in ARMS {
            for position in 0..4u64 {
                let id = format!("r{round}_{arm}_t{position}");
                let exit = fs::read_to_string(root.join(format!("{id}.exit"))).unwrap();
                assert_eq!(exit.trim(), "0", "{id}");
                let row: Value = serde_json::from_slice(
                    &fs::read(root.join(format!("{id}.json"))).unwrap()).unwrap();
                let meter: Value = serde_json::from_slice(
                    &fs::read(root.join(format!("{id}.meter.json"))).unwrap()).unwrap();
                assert_eq!(row["schema"], if ["inverse-velu", "full-image"].contains(&arm) {
                    "different_j_chi_heldout_run/v4"
                } else {
                    "different_j_chi_heldout_run/v3"
                });
                assert_eq!(meter["schema"], "different_j_chi_heldout_meter/v1");
                assert_eq!(row["verified"], true);
                assert_eq!(meter["exit_code"], 0);
                assert_eq!(meter["timed_out"], false);
                assert_eq!(row["round"], round);
                assert_eq!(row["position"], position);
                assert_eq!(row["arm"], arm);
                let g = group_count(&row);
                let cpu_us = number(&meter, "child_user_cpu_us")
                    + number(&meter, "child_system_cpu_us");
                let wall_ns = number(&meter, "wall_ns");
                let this = totals.entry(arm.to_string()).or_default();
                this.count += 1;
                this.group += g;
                this.meter_cpu_us += cpu_us;
                this.meter_wall_ns += wall_ns;
                this.solver_cpu_ns += number(&row, "total_cpu_ns");
                this.peak_rss_bytes = this.peak_rss_bytes.max(number(&meter, "child_peak_rss_bytes"));
                this.trials += number(&row, "trials");
                this.relations += number(&row, "relations");
                this.rank += number(&row, "rank");
                this.fp_muls += ["setup_fp_muls", "query_fp_muls", "verify_fp_muls"]
                    .iter().map(|key| number(&row, key)).sum::<u64>();
                this.fp_invs += ["setup_fp_invs", "query_fp_invs", "verify_fp_invs"]
                    .iter().map(|key| number(&row, key)).sum::<u64>();
                paired.entry((round, position)).or_default()
                    .insert(arm.to_string(), (g, cpu_us, number(&row, "trials"),
                        number(&row, "relations")));
            }
        }
    }
    assert_eq!(paired.len(), 12);
    for arms in paired.values() {
        assert_eq!(arms.len(), 9);
        let control = arms["direct-pair"];
        for arm in ["chi-resultant", "source-s3", "quotient-s3",
            "inverse-velu", "full-image"] {
            assert_eq!(arms[arm].2, control.2, "trial support differs in {arm}");
            assert_eq!(arms[arm].3, control.3, "relation support differs in {arm}");
        }
    }
    let mut rows = Vec::new();
    for arm in ARMS {
        let t = &totals[arm];
        assert_eq!(t.count, 12);
        rows.push(json!({
            "arm": arm,
            "children": t.count,
            "mean_cold_group_s": mean(t.group,t.count)/(229f64).sqrt(),
            "mean_cold_group_ops": mean(t.group,t.count),
            "mean_whole_child_cpu_ms": mean(t.meter_cpu_us,t.count)/1000.0,
            "mean_whole_child_wall_ms_unisolated": mean(t.meter_wall_ns,t.count)/1_000_000.0,
            "mean_solver_thread_cpu_ms_diagnostic": mean(t.solver_cpu_ns,t.count)/1_000_000.0,
            "max_peak_rss_bytes": t.peak_rss_bytes,
            "mean_trials": mean(t.trials,t.count),
            "mean_relations": mean(t.relations,t.count),
            "mean_rank": mean(t.rank,t.count),
            "mean_reported_fp_muls_overlaps_group": mean(t.fp_muls,t.count),
            "mean_reported_fp_invs_overlaps_group": mean(t.fp_invs,t.count)
        }));
    }
    let rho = &totals["rho-negation"];
    let direct = &totals["direct-pair"];
    let ratios: Vec<_> = ["chi-resultant", "direct-pair", "source-s3", "quotient-s3",
        "inverse-velu", "full-image",
        "source-bsgs", "rho-negation-aa"].iter().map(|&arm| {
        let t = &totals[arm];
        json!({"arm":arm,
            "mean_group_vs_rho":t.group as f64/rho.group as f64,
            "mean_whole_child_cpu_vs_rho":t.meter_cpu_us as f64/rho.meter_cpu_us as f64,
            "mean_group_vs_direct_pair":t.group as f64/direct.group as f64,
            "mean_whole_child_cpu_vs_direct_pair":t.meter_cpu_us as f64/direct.meter_cpu_us as f64})
    }).collect();
    println!("{}",json!({
        "schema":"different_j_chi_heldout_summary/v4",
        "status":"pass",
        "p":211,"r":229,"paired_round_targets":12,"children":108,
        "rows":rows,"ratios":ratios,
        "cost_scope":"cold whole-child group counts and metered CPU; field counters overlap group; wall unisolated; portable S remains unset"
    }));
}
