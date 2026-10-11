//! First distributional relation-yield probe for Koblitz (binary) fixtures.
//!
//! Freezes n=19, a=1, signed_expanded, and sweeps η numerators against a
//! fixed denominator, reporting exact-coverage hit counts from the producer.
//!
//! ```bash
//! cargo build --release --example koblitz_rank_fixture
//! cargo run --release --example relation_yield_koblitz_eta -- 19 1 8
//! ```
//!
//! Public synthetic only.

use std::process::Command;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(19);
    let a: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let eta_den: u32 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(8);
    let bin = std::env::var("KIC_RANK_BIN").unwrap_or_else(|_| {
        "target/release/examples/koblitz_rank_fixture".into()
    });

    let mut rows = Vec::new();
    for eta_num in 1..=eta_den {
        let t0 = Instant::now();
        let out = Command::new(&bin)
            .args([
                &n.to_string(),
                &a.to_string(),
                &eta_num.to_string(),
                &eta_den.to_string(),
                "1",
                "signed_expanded",
                "independent",
                "pair_pair_16",
                "1",
            ])
            .env("KIC_SUMMARY_ONLY", "1")
            .env("KIC_EXACT_COVERAGE_ONLY", "1")
            .env("KIC_INCREMENTAL_RANK_CROSSCHECK", "0")
            .output()
            .expect("spawn koblitz_rank_fixture");
        let wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        if !out.status.success() {
            rows.push(serde_json::json!({
                "eta_num": eta_num,
                "eta_den": eta_den,
                "ok": false,
                "stderr": String::from_utf8_lossy(&out.stderr).chars().take(400).collect::<String>(),
            }));
            continue;
        }
        let text = String::from_utf8_lossy(&out.stdout);
        let obj: serde_json::Value = text
            .lines()
            .find_map(|l| serde_json::from_str(l).ok())
            .or_else(|| serde_json::from_str(text.trim()).ok())
            .unwrap_or(serde_json::json!({}));
        let reach = obj.get("reachable_nonzero_subgroup_targets").and_then(|v| v.as_u64());
        let complete = obj.get("complete_subgroup_coverage").and_then(|v| v.as_bool());
        let f = obj.get("factor_base_points").and_then(|v| v.as_u64());
        let pairs = obj.get("distinct_pair_sum_points").and_then(|v| v.as_u64());
        let base_hash = obj.get("base_hash").and_then(|v| v.as_str()).unwrap_or("");
        // Hit-rate proxy: reachable / (2^n - 1) when present; else null.
        let subgroup_nontrivial = (1u64 << n).saturating_sub(1);
        let hit = reach.map(|r| r as f64 / subgroup_nontrivial as f64);
        rows.push(serde_json::json!({
            "eta_num": eta_num,
            "eta_den": eta_den,
            "ok": true,
            "factor_base_points": f,
            "reachable_nonzero_subgroup_targets": reach,
            "complete_subgroup_coverage": complete,
            "distinct_pair_sum_points": pairs,
            "hit_rate_proxy": hit,
            "base_hash": base_hash,
            "wall_ms": wall_ms,
        }));
        eprintln!(
            "eta={eta_num}/{eta_den} reach={reach:?} complete={complete:?} wall_ms={wall_ms:.1}"
        );
    }

    let report = serde_json::json!({
        "schema_version": 2,
        "stage": "relation_yield",
        "regime": "binary",
        "n_or_bits": n,
        "a": a,
        "base_id_or_hash": "koblitz_rank_fixture signed_expanded independent pair_pair_16 seed=1",
        "eta_or_coverage_policy": format!("sweep eta_num=1..{eta_den} / {eta_den}"),
        "pr_decomposition_or_hit_rate_with_ci": {
            "definition": "reachable_nonzero_subgroup_targets / (2^n - 1)",
            "ci": "not_computed_single_draw_per_eta",
            "rows": rows,
        },
        "trials_per_relation": "not_applicable_exact_coverage_scan",
        "target_mix": {
            "natural": "all_nonzero_subgroup_via_exact_coverage",
            "planted_sat": 0,
            "proven_unsat": 0,
        },
        "claim_boundary": "Public synthetic Koblitz exact-coverage eta sweep only.",
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic yield law",
            "not ledger promotion until CI/distributional gates",
        ],
    });
    println!("{}", serde_json::to_string_pretty(&report).unwrap());
}
