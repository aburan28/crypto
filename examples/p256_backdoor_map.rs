//! Print the P-256 backdoor map as Markdown and re-run the X9.62 seed →
//! `b` derivation against the published constants.
//!
//! ```bash
//! cargo run --release --example p256_backdoor_map            # Markdown to stdout
//! cargo run --release --example p256_backdoor_map -- --json  # machine-readable summary
//! ```
//!
//! Exits non-zero if the derivation check fails, so the example doubles
//! as a regression guard for the seed facts the note relies on.

use crypto_lib::cryptanalysis::p256_backdoor_map::{
    format_backdoor_assessment, p256_backdoor_assessment, verify_p256_seed_derivation,
    GrindCostModel, P256_X962_R_HEX,
};

fn main() {
    let json = std::env::args().any(|a| a == "--json");
    let derivation = verify_p256_seed_derivation();
    let assessment = p256_backdoor_assessment();
    if json {
        let hyps: Vec<serde_json::Value> = assessment
            .hypotheses
            .iter()
            .map(|h| {
                serde_json::json!({
                    "rank": h.rank,
                    "locus": format!("{:?}", h.locus),
                    "name": h.name,
                    "grind_model": format!("{:?}", h.grind),
                    "log2_candidates_1999": assessment.budget.log2_candidates(h.grind),
                    "public_status": format!("{:?}", h.public_status),
                    "verdict": format!("{:?}", h.verdict),
                    "repo_probe": h.repo_probe,
                })
            })
            .collect();
        let doc = serde_json::json!({
            "status": format!("{:?}", assessment.status),
            "seed_hex": assessment.seed.seed_hex,
            "x962": {"s": derivation.s, "v": derivation.v, "r_hex": format!("{:x}", derivation.r),
                      "r_pinned_hex": P256_X962_R_HEX, "holds": derivation.holds},
            "dual_ec_implicates_curve_math": assessment.confirmed.implicates_curve_math,
            "grind_budget": {"cpu_years": assessment.budget.cpu_years,
                "log2_coefficient_tests": assessment.budget.log2_candidates(GrindCostModel::CoefficientTest),
                "log2_point_counts": assessment.budget.log2_candidates(GrindCostModel::PointCount)},
            "hypotheses": hyps,
        });
        println!("{}", serde_json::to_string_pretty(&doc).unwrap());
    } else {
        print!("{}", format_backdoor_assessment(&assessment));
    }
    if !derivation.holds || format!("{:x}", derivation.r) != P256_X962_R_HEX {
        eprintln!("X9.62 derivation check FAILED");
        std::process::exit(1);
    }
}
