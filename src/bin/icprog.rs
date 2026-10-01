//! `icprog`: the `ic` tool programme's harness, native.
//!
//! AGENTS.md asks that research tooling be native: experiment harnesses,
//! benchmark drivers, correctness tests, replay tools and the generation
//! of research results.  The programme's harness began in Python
//! (`research/ic_tool_program/harness/`, each round's `run.py` and
//! `analyse.py`); this binary replaces it piece by piece, each piece
//! checked against the frozen outputs of the one it replaces.
//!
//! It is a binary of its own, apart from `ic`, so that the harness never
//! changes the code it measures.
//!
//! - `icprog analyse <round>`: a round's figures and decision, from its run
//!   tree only.  `r05` is the native analysis of R05; `r03` reproduces
//!   R03's frozen `analysis.json` from its frozen runs.
// Shared with `isolated_bench`, which uses parts this binary does not.
#[path = "icprog/json.rs"]
#[allow(dead_code)]
mod json;
#[path = "icprog/rounds.rs"]
mod rounds;
#[path = "icprog/runs.rs"]
mod runs;
#[path = "icprog/stats.rs"]
mod stats;
#[path = "icprog/suite.rs"]
mod suite;

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, Subcommand, ValueEnum};

#[derive(Parser)]
#[command(
    name = "icprog",
    version,
    about = "The ic tool programme's harness: round analyses, natively"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, ValueEnum)]
enum Round {
    R03,
    R05,
}

#[derive(Subcommand)]
enum Command {
    /// A round's figures and decision, from its run tree only, as JSON.
    Analyse {
        round: Round,
        /// The repository checkout (default: the current directory).
        #[arg(long, default_value = ".")]
        root: PathBuf,
        /// The run tree (default: the round's `runs/`).
        #[arg(long)]
        runs: Option<PathBuf>,
    },
}

fn analyse(round: Round, root: PathBuf, runs: Option<PathBuf>) -> Result<String, String> {
    let programme = suite::programme(&root)?;
    let dir = match round {
        Round::R03 => "R03-curve-construction",
        Round::R05 => "R05-presence-filter",
    };
    let round_dir = programme.join("rounds").join(dir);
    let runs = runs.unwrap_or_else(|| round_dir.join("runs"));
    if !runs.is_dir() {
        return Err(format!("{} is not a run tree", runs.display()));
    }
    let ctx = rounds::Ctx {
        programme,
        round_dir,
        runs,
    };
    let doc = match round {
        Round::R03 => rounds::r03::analyse(&ctx)?,
        Round::R05 => rounds::r05::analyse(&ctx)?,
    };
    Ok(json::dumps(&doc, 1))
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Analyse { round, root, runs } => analyse(round, root, runs),
    };
    match result {
        Ok(text) => {
            println!("{text}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("icprog: {e}");
            ExitCode::FAILURE
        }
    }
}
