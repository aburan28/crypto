//! Run the round-3 P-256 trapdoor probes: the honest-distribution
//! comparison of the published constants against simulated X9.62
//! candidates, plus the seed-preimage phrase probe.
//!
//! ```bash
//! cargo run --release --example p256_trapdoor_round3 -- 1048576
//! ```

use crypto_lib::cryptanalysis::p256_trapdoor_round3::{
    curated_phrase_list, p256_trapdoor_round3_report, seed_preimage_probe, P256_SEED_HEX,
};

fn main() {
    let n_seeds: u64 = std::env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(1024);

    eprintln!("simulating {n_seeds} X9.62 candidates (deterministic splitmix seeds)...");
    let records = p256_trapdoor_round3_report(n_seeds);

    println!("# P-256 trapdoor round 3: honest-distribution comparison\n");
    println!("seed count: {n_seeds}; probes: Bonferroni-corrected\n");
    for r in &records {
        println!("- **{}**", r.name);
        println!("  - published: {}", r.published);
        println!("  - honest:   {}", r.honest);
        println!("  - p (corrected): {:.6}", r.p_value);
    }

    let phrases = curated_phrase_list();
    let preimage = seed_preimage_probe(&phrases);
    println!(
        "\n# Seed-preimage phrase probe\n\n- target: {P256_SEED_HEX}\n- phrases tried: {}\n- result: {}",
        phrases.len(),
        preimage.as_deref().unwrap_or("no preimage found")
    );

    if records.iter().any(|r| r.p_value < 0.001) {
        eprintln!("round-3 outlier flagged");
        std::process::exit(1);
    }
}
