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
//! - `icprog run <round> <step>`: a round's declared timed steps on the
//!   native runner (`harness/bench.py`, ported), through `isolated_bench`.
//!   It resumes a run tree where it stopped and never overwrites.
#[path = "icprog/bench.rs"]
mod bench;
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
    /// One of a round's declared timed steps, natively; R05 for now.
    Run {
        round: Round,
        /// plan, manifest-resumed, compare, holdout or extend.
        step: String,
        /// The repository checkout (default: the current directory).
        #[arg(long, default_value = ".")]
        root: PathBuf,
        /// The run tree (default: the round's `runs/`).
        #[arg(long)]
        runs: Option<PathBuf>,
        /// The baseline arm's `ic`.
        #[arg(long)]
        base: PathBuf,
        /// The candidate arm's `ic`.
        #[arg(long)]
        cand: PathBuf,
        /// The isolation tool (default: `isolated_bench` beside this binary).
        #[arg(long)]
        isolate: Option<PathBuf>,
    },
}

fn context(
    round: Round,
    root: &std::path::Path,
    runs: Option<PathBuf>,
) -> Result<rounds::Ctx, String> {
    let programme = suite::programme(root)?;
    let dir = match round {
        Round::R03 => "R03-curve-construction",
        Round::R05 => "R05-presence-filter",
    };
    let round_dir = programme.join("rounds").join(dir);
    let runs = runs.unwrap_or_else(|| round_dir.join("runs"));
    if !runs.is_dir() {
        return Err(format!("{} is not a run tree", runs.display()));
    }
    Ok(rounds::Ctx {
        programme,
        round_dir,
        runs,
    })
}

fn run(
    round: Round,
    step: &str,
    root: PathBuf,
    runs: Option<PathBuf>,
    arms: [PathBuf; 2],
    isolate: Option<PathBuf>,
) -> Result<String, String> {
    let ctx = context(round, &root, runs)?;
    let isolate = match isolate {
        Some(p) => p,
        None => std::env::current_exe()
            .map_err(|e| e.to_string())?
            .with_file_name("isolated_bench"),
    };
    if !isolate.exists() {
        return Err(format!("no isolation tool at {}", isolate.display()));
    }
    let [base, cand] = arms;
    for b in [&base, &cand] {
        if !b.exists() {
            return Err(format!("no binary at {}", b.display()));
        }
    }
    let arms = [
        bench::Arm {
            name: "base".into(),
            binary: base,
        },
        bench::Arm {
            name: "cand".into(),
            binary: cand,
        },
    ];
    let b = bench::Bench { isolate };
    let doc = match round {
        Round::R05 => rounds::r05::run(&ctx, step, &b, &arms, &root)?,
        Round::R03 => return Err("R03 is complete; its runs are frozen".into()),
    };
    Ok(match doc {
        json::J::Null => String::new(),
        d => json::dumps(&d, 1),
    })
}

fn analyse(round: Round, root: PathBuf, runs: Option<PathBuf>) -> Result<String, String> {
    let ctx = context(round, &root, runs)?;
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
        Command::Run {
            round,
            step,
            root,
            runs,
            base,
            cand,
            isolate,
        } => run(round, &step, root, runs, [base, cand], isolate),
    };
    match result {
        Ok(text) if text.is_empty() => ExitCode::SUCCESS,
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
