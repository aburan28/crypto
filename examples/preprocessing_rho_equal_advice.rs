//! Lead A of `RESEARCH_IC_LEADS_AND_CONNECTIONS_20261010.md` §2, as a
//! calculation from recorded values: the compact-orbit IC "online after
//! precompute" target stage against Bernstein–Lange preprocessing rho with
//! the same advice (stored entries), on the n = 73 and n = 83 single-target
//! fixtures.
//!
//! Model (the note's): folded group size `N = r / (2n)`; advice `S` = the IC
//! root-table entries; BL online `T ≈ 2·√(N/S)` walk steps, BL precompute
//! `≈ √(N·S)` steps; the generic preprocessing line `S·T² ≳ N`. Rho walk
//! rates (steps/s) are the recorded paired rho walks on the same host, so
//! the BL time is a projection, not a measurement.
//!
//! Usage: `preprocessing_rho_equal_advice <inputs.json>` (prints a Markdown
//! table; `research/ic_preprocessing_rho_equal_advice_20261011/inputs.json`).
use serde_json::Value;

fn f(v: &Value) -> f64 {
    v.as_f64().expect("number")
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: <inputs.json>");
    let d: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    println!("| rung | log₂ N | log₂ S | BL online steps | IC target probes | ops ratio IC/BL | rho rate (steps/s, R1–R3) | BL online ms (projected) | IC target ms (measured) | log₂ S·T² (IC) − log₂ N | BL precompute steps | IC rank-stage probes (600 × mean R1) | plain rho ideal steps |");
    println!("|:--|--:|--:|--:|--:|--:|:--|:--|:--|--:|--:|--:|--:|");
    for key in ["n73", "n83"] {
        let v = &d[key];
        let n: f64 = key[1..].parse().unwrap();
        let log_n_fold = f(&v["log2_r"]) - (2.0 * n).log2();
        let s = f(&v["root_table_entries"]);
        let log_s = s.log2();
        let t_bl = 2.0 * ((log_n_fold - log_s) / 2.0).exp2();
        let probes = f(&v["ic_target_probes"][0]);
        let rates: Vec<f64> = (0..3)
            .map(|i| f(&v["rho_walk_steps"][i]) / (f(&v["rho_walk_ms"][i]) / 1000.0))
            .collect();
        let bl_ms: Vec<String> = rates
            .iter()
            .map(|r| format!("{:.1}", t_bl / r * 1000.0))
            .collect();
        let ic_ms: Vec<String> = (0..3)
            .map(|i| format!("{:.0}", f(&v["ic_target_ms"][i])))
            .collect();
        let above = log_s + 2.0 * probes.log2() - log_n_fold;
        let bl_pre = ((log_n_fold + log_s) / 2.0).exp2();
        let ic_pre = 600.0 * f(&v["ic_rank_probes_mean"][0]);
        println!(
            "| {key} | {log_n_fold:.2} | {log_s:.2} | {t_bl:.0} | {probes:.0} | {:.0} | {} | {} | {} | {above:.1} | {bl_pre:.2e} | {ic_pre:.2e} | {:.2e} |",
            probes / t_bl,
            rates.iter().map(|r| format!("{r:.0}")).collect::<Vec<_>>().join(", "),
            bl_ms.join(", "),
            ic_ms.join(", "),
            f(&v["rho_ideal_steps"]),
        );
    }
}
