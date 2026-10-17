//! A small read-only consumer of the ecbench evidence format. It intentionally
//! builds independently of the cryptography library. Pricing and bootstrap
//! implementations are shared with the harness, not independently reimplemented.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;

mod evidence;
mod leaderboard;

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
}

fn run(cli: Cli) -> Result<(), String> {
    let board = leaderboard::build(&cli.root, &cli.manifest)?;
    leaderboard::write(&board, &cli.out, cli.check)?;
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
