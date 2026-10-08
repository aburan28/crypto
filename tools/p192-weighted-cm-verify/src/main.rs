use std::{path::PathBuf, process::ExitCode};

use clap::{Parser, Subcommand};
use p192_weighted_cm_verify::verify::run_preflight;

#[derive(Debug, Parser)]
#[command(name = "p192_weighted_cm_verify")]
#[command(about = "Independent P-192 weighted-CM artifact verifier")]
struct Arguments {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Preflight {
        #[arg(long)]
        source: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
}

fn main() -> ExitCode {
    let arguments = Arguments::parse();
    let result = match arguments.command {
        Command::Preflight { source, out } => run_preflight(&source, &out),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("p192_weighted_cm_verify: {error}");
            ExitCode::FAILURE
        }
    }
}
