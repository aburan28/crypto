//! Native P-256 affine-bitbox factor-base inventory.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_bitbox_factor_base;

#[derive(Parser)]
#[command(about = "Build the frozen registered P-256 affine-bitbox factor base")]
struct Cli {
    /// Registered ICV1 slug.  The exact P-256 model is the only accepted curve.
    #[arg(long)]
    curve: String,
    /// Compact result JSON path; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Optional full `ecbench.factor_base_dump/v1-wide` JSON path.
    #[arg(long)]
    dump_out: Option<PathBuf>,
    /// Rebuild all eight candidates and compare every emitted point row.
    #[arg(long)]
    verify: bool,
}

fn write_json(path: Option<&PathBuf>, value: &impl serde::Serialize) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|error| error.to_string())? + "\n";
    match path {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string()),
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

fn run(cli: Cli) -> Result<(), String> {
    if cli.curve != crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG {
        return Err(format!(
            "P-256 affine-bitbox builder has no registered curve `{}`",
            cli.curve
        ));
    }
    let result = p256_bitbox_factor_base::build()?;
    if cli.verify {
        p256_bitbox_factor_base::verify(&result)?;
    }
    if let Some(path) = cli.dump_out.as_ref() {
        write_json(Some(path), &result.build.dump)?;
    }
    let compact = p256_bitbox_factor_base::compact(&result, cli.verify);
    write_json(cli.out.as_ref(), &compact)?;
    eprintln!(
        "{} on {}: {} points, {} columns, offset {}{}",
        compact.fb_id,
        compact.curve,
        compact.signed_points,
        compact.columns,
        compact.selected_offset,
        if compact.verified { " (verified)" } else { "" },
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
