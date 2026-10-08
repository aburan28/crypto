//! Native runner and verifier for `EXP-SCURVE-29040c`.

use clap::{Args, Parser, Subcommand};
use crypto_lib::cryptanalysis::p192_singular_recovery::{
    run_recovery, verify_report, RecoveryConfig, RecoveryReport,
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(about = "P-192 singular-companion x-only planted-key recovery")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Execute the frozen native recovery and write its JSON certificate.
    Run(RunArgs),
    /// Replay and verify a previously written JSON certificate.
    Verify {
        #[arg(long)]
        input: PathBuf,
    },
}

#[derive(Args)]
struct RunArgs {
    #[arg(long)]
    seed_hex: String,
    #[arg(long)]
    oracle: String,
    #[arg(long)]
    residue_mode: String,
    #[arg(long)]
    orientation_search: String,
    #[arg(long, default_value_t = 4)]
    threads: usize,
    #[arg(long, default_value_t = 256)]
    batch_lanes: usize,
    #[arg(long, default_value_t = 24)]
    table_bits: u32,
    #[arg(long)]
    output: PathBuf,
}

fn write_report(path: &Path, report: &RecoveryReport) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("failed to create {}: {error}", parent.display()))?;
    }
    let bytes = serde_json::to_vec_pretty(report)
        .map_err(|error| format!("failed to serialize recovery report: {error}"))?;
    fs::write(path, bytes).map_err(|error| format!("failed to write {}: {error}", path.display()))
}

fn run(args: RunArgs) -> Result<(), String> {
    if args.oracle != "x-only-hmac-sha256" {
        return Err("--oracle must be x-only-hmac-sha256".into());
    }
    if args.residue_mode != "nested-composite" {
        return Err("--residue-mode must be nested-composite".into());
    }
    if args.orientation_search != "interleaved" {
        return Err("--orientation-search must be interleaved".into());
    }
    let seed = hex::decode(&args.seed_hex)
        .map_err(|error| format!("--seed-hex is not valid hexadecimal: {error}"))?;
    let config = RecoveryConfig::frozen(seed, args.threads, args.batch_lanes, args.table_bits);
    let report = run_recovery(&config)?;
    write_report(&args.output, &report)?;
    println!(
        "recovered d={} certificate_sha256={} output={}",
        report.interval.recovered_d,
        report.certificate_sha256,
        args.output.display()
    );
    Ok(())
}

fn verify(input: &Path) -> Result<(), String> {
    let bytes =
        fs::read(input).map_err(|error| format!("failed to read {}: {error}", input.display()))?;
    let report: RecoveryReport = serde_json::from_slice(&bytes)
        .map_err(|error| format!("failed to parse {}: {error}", input.display()))?;
    verify_report(&report)?;
    println!(
        "verified certificate_sha256={} recovered_d={}",
        report.certificate_sha256, report.interval.recovered_d
    );
    Ok(())
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Run(args) => run(args),
        Command::Verify { input } => verify(&input),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p192_singular_recovery: {error}");
            ExitCode::FAILURE
        }
    }
}
