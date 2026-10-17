//! A small read-only consumer of the ecbench evidence format. It intentionally
//! builds independently of the cryptography library. Pricing and bootstrap
//! implementations are shared with the harness, not independently reimplemented.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand};

mod evidence;
mod leaderboard;
mod reproduction;

#[allow(dead_code)]
#[path = "../../../src/cryptanalysis/ecbench/band.rs"]
mod band;
#[allow(dead_code)]
#[path = "../../../src/cryptanalysis/ecbench/canonical.rs"]
mod canonical;
#[allow(dead_code)]
#[path = "../../../src/cryptanalysis/ecbench/stats.rs"]
mod stats;

// Minimal API adapters for the reused, read-only source modules above.
mod methods {
    include!(concat!(env!("OUT_DIR"), "/priced_native.rs"));
}
mod record {
    pub use crate::evidence::Record;
}
mod hash {
    pub mod sha256 {
        pub fn sha256(bytes: &[u8]) -> [u8; 32] {
            use sha2::Digest;
            sha2::Sha256::digest(bytes).into()
        }
    }
}
mod cryptanalysis {
    pub mod ecbench {
        pub mod canonical {
            pub use crate::canonical::splitmix64;
        }
    }
    pub mod ic_boundary {
        pub use crate::evidence::pinned_calibration_ratios;
    }
}

#[derive(Parser)]
#[command(
    about = "Backend coverage and observed frontiers from sealed ecbench sessions; no solver execution"
)]
struct Cli {
    #[arg(long, default_value = "docs/backend-leaderboard/manifest.json")]
    manifest: PathBuf,
    #[arg(long, default_value = ".")]
    root: PathBuf,
    #[arg(long, default_value = "docs/backend-leaderboard")]
    out: PathBuf,
    /// Check all generated views without writing them.
    #[arg(long)]
    check: bool,
    /// Explicitly replace the reviewed input/output lock after a deliberate update.
    #[arg(long, conflicts_with = "check")]
    seal: bool,
    #[arg(
        long,
        default_value = "docs/backend-leaderboard/reproduction.lock.json"
    )]
    lock: PathBuf,
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand)]
enum Command {
    /// Recompute every view and verify its pinned hash; executes no solver.
    Reproduce {
        /// Obtain this hash from a trusted, reviewed commit, not this working tree.
        #[arg(long)]
        expected_lock_sha256: String,
    },
    /// Compare a freshly audited session to one exact frozen source, including solver statistics.
    CheckRerun {
        /// Session directory exactly as named in the manifest.
        #[arg(long)]
        source: PathBuf,
        #[arg(long)]
        rerun: PathBuf,
        /// Fresh --replay-all audit of the frozen original.
        #[arg(long)]
        original_audit: PathBuf,
        /// Fresh --replay-all audit of the new session.
        #[arg(long)]
        rerun_audit: PathBuf,
        /// The actual ecbench executable used for the run and both audits.
        #[arg(long)]
        binary: PathBuf,
        #[arg(long)]
        expected_lock_sha256: String,
    },
}

fn run(cli: Cli) -> Result<(), String> {
    if (cli.check || cli.seal) && cli.command.is_some() {
        return Err("--check/--seal cannot be combined with a subcommand".into());
    }
    let lock = if cli.seal {
        reproduction::compiled_check(&cli.root)?;
        None
    } else {
        let expected = match &cli.command {
            Some(Command::Reproduce {
                expected_lock_sha256,
            })
            | Some(Command::CheckRerun {
                expected_lock_sha256,
                ..
            }) => Some(expected_lock_sha256.as_str()),
            None => None,
        };
        Some(reproduction::verify_inputs(
            &cli.root,
            &cli.manifest,
            &cli.lock,
            expected,
        )?)
    };
    let board = leaderboard::build(&cli.root, &cli.root.join(&cli.manifest))?;
    if let Some(lock) = &lock {
        reproduction::verify_outputs(lock, &board)?;
    }
    if let Some(command) = cli.command {
        let mut receipt = match command {
            Command::Reproduce { .. } => {
                reproduction::report_receipt(&cli.lock, lock.as_ref().unwrap(), &board)?
            }
            Command::CheckRerun {
                source,
                rerun,
                original_audit,
                rerun_audit,
                binary,
                ..
            } => reproduction::check_rerun(
                &cli.root,
                &board,
                &source,
                &rerun,
                &original_audit,
                &rerun_audit,
                &binary,
            )?,
        };
        receipt["lock_sha256"] = serde_json::Value::String(canonical::sha256_hex(
            &std::fs::read(&cli.lock).map_err(|e| e.to_string())?,
        ));
        println!(
            "{}",
            serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())?
        );
        return Ok(());
    }
    leaderboard::write(&board, &cli.out, cli.check)?;
    if cli.seal {
        reproduction::seal(&cli.root, &cli.manifest, &cli.lock, &board)?;
    }
    eprintln!("{} frozen sessions, {} exact-target backend entries, {} planned cells without measurements; {}",
        board.sources.len(),board.entries.len(),board.unattempted.len(),if cli.check {"current"} else {"written"});
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
