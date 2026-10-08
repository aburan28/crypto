//! Native wide-integer factor-base inventory for the registered P-256 model.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base;

#[derive(Parser)]
#[command(about = "Build and verify the registered P-256 Dickson factor base")]
struct Cli {
    /// Registered ICV1 slug.  The exact P-256 model is the only accepted curve.
    #[arg(long)]
    curve: String,
    /// Factor-base specification, for example `dickson-torus:depth=18`.
    #[arg(long)]
    factor_base: String,
    /// JSON `ecbench.factor_base_dump/v1-wide` output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Optional SQL loader for the repository's canonical factor-base tables.
    #[arg(long)]
    sql_out: Option<PathBuf>,
    /// Report the exact signed domain and Poisson heuristic at this length.
    #[arg(long)]
    relation_length: Option<u32>,
    /// Rebuild and compare every emitted identity and point row.
    #[arg(long)]
    verify: bool,
}

fn run(cli: Cli) -> Result<(), String> {
    if cli.curve != p256_dickson_factor_base::CURVE_SLUG {
        return Err(format!(
            "P-256 factor-base builder has no registered curve `{}`; known: {}",
            cli.curve,
            p256_dickson_factor_base::CURVE_SLUG
        ));
    }
    let built = p256_dickson_factor_base::build(&cli.factor_base)?;
    if cli.verify {
        p256_dickson_factor_base::verify(&built.dump)?;
    }
    if let Some(length) = cli.relation_length {
        let metrics =
            p256_dickson_factor_base::relation_metrics(built.dump.factor_base.columns, length)?;
        eprintln!(
            "relation m={}: signed_domain={}, lambda={:.15}, success={:.12}%, selector ideal regularity <= {}",
            metrics.relation_length,
            metrics.signed_domain,
            metrics.poisson_mean,
            100.0 * metrics.poisson_success,
            metrics.selector_ideal_regularity_bound,
        );
    }
    let text = serde_json::to_string_pretty(&built.dump).map_err(|e| e.to_string())? + "\n";
    match &cli.out {
        Some(path) => std::fs::write(path, &text).map_err(|e| e.to_string())?,
        None => print!("{text}"),
    }
    if let Some(path) = &cli.sql_out {
        let sql = p256_dickson_factor_base::sql_for_dump(&built.dump)?;
        std::fs::write(path, sql).map_err(|e| e.to_string())?;
    }
    eprintln!(
        "{} on {}: {} points, {} columns, legacy digest {}{}{}",
        built.dump.factor_base.fb_id,
        built.dump.curve.slug,
        built.dump.factor_base.signed_points,
        built.dump.factor_base.columns,
        built.legacy_enumeration_sha256,
        cli.out
            .as_ref()
            .map(|path| format!(" -> {}", path.display()))
            .unwrap_or_default(),
        if cli.verify { " (verified)" } else { "" },
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
