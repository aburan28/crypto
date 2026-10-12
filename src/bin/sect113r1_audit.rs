//! Native producer and replay verifier for the exact sect113r1 audit.

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::sect113r1_audit::{run_audit, verify_report, AuditReport};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(about = "Exact sect113r1 class, isogeny, and invalid-input audit")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Execute the deterministic audit and write its typed JSON certificate.
    Run {
        #[arg(long)]
        output: PathBuf,
    },
    /// Recompute every certificate and compare it with a saved report.
    Verify {
        #[arg(long)]
        input: PathBuf,
    },
}

fn write_report(path: &Path, report: &AuditReport) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("failed to create {}: {error}", parent.display()))?;
    }
    let bytes = serde_json::to_vec_pretty(report)
        .map_err(|error| format!("failed to serialize report: {error}"))?;
    fs::write(path, bytes).map_err(|error| format!("failed to write {}: {error}", path.display()))
}

fn run(output: &Path) -> Result<(), String> {
    let report = run_audit()?;
    write_report(output, &report)?;
    println!(
        "diagnostic={} class={} implementation={} isogenies={} sha256={} output={}",
        report.evidence_status,
        report.class_weakness.verdict,
        report.singular_companion.verdict,
        report.valid_isogenies.rational_degree_five_kernel_count,
        report.certificate_sha256,
        output.display()
    );
    Ok(())
}

fn verify(input: &Path) -> Result<(), String> {
    let bytes =
        fs::read(input).map_err(|error| format!("failed to read {}: {error}", input.display()))?;
    let report: AuditReport = serde_json::from_slice(&bytes)
        .map_err(|error| format!("failed to parse {}: {error}", input.display()))?;
    verify_report(&report)?;
    println!(
        "verified diagnostic={} sha256={}",
        report.evidence_status, report.certificate_sha256
    );
    Ok(())
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Run { output } => run(&output),
        Command::Verify { input } => verify(&input),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("sect113r1_audit: {error}");
            ExitCode::FAILURE
        }
    }
}
