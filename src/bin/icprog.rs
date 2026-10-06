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
//!   tree only.  `r05`, `r02b` and `r07` are native; `r03` reproduces R03's
//!   frozen `analysis.json` from its frozen runs.
//! - `icprog holdouts`: a round's fresh holdouts by suite v1's own
//!   construction (`suite/v1/make_suite.py`, ported), drawn once, or
//!   checked against the files a round froze.
//! - `icprog run <round> <step>`: a round's declared steps on the native
//!   runner (`harness/bench.py`, ported), through `isolated_bench`: its
//!   manifest, pin, A/A, timed rows and callgrind profiles.  It resumes a
//!   run tree where it stopped and never overwrites.
//! - `icprog table <round> <analysis.json>`: the README's tables, rendered
//!   from the analysis only.
//! - `icprog baseline`: an accepted round's ledger entry, appended.
//! - `icprog callgrind-phases` and `icprog callgrind-control`: a callgrind
//!   profile split by phase (`harness/callgrind_phases.py`, ported), and
//!   R02b's control over a run tree's profiles.
//! - `icprog rule <step>`: the single-target rule's comparison (ledger
//!   §23's `run.py`, `claims.py` and `analyse.py`, ported): its rows on the
//!   native runner, each row's identities, independent replay and checked
//!   claim, and the figures read from the checked rows.
#[path = "icprog/autolab.rs"]
mod autolab;
#[path = "icprog/bench.rs"]
mod bench;
#[path = "icprog/callgrind.rs"]
mod callgrind;
#[path = "icprog/f5_control.rs"]
mod f5_control;
#[path = "icprog/f5_control_publication.rs"]
mod f5_control_publication;
#[path = "icprog/f5_target.rs"]
mod f5_target;
#[path = "icprog/identity.rs"]
mod identity;
#[path = "icprog/oracle.rs"]
mod oracle;
#[path = "icprog/report.rs"]
mod report;
// Shared with `isolated_bench`, which uses parts this binary does not.
#[path = "icprog/json.rs"]
#[allow(dead_code)]
mod json;
#[path = "icprog/pin.rs"]
mod pin;
#[path = "icprog/prepared_f5.rs"]
mod prepared_f5;
#[path = "icprog/prepared_sat.rs"]
mod prepared_sat;
#[path = "icprog/pyrandom.rs"]
mod pyrandom;
#[path = "icprog/rounds.rs"]
mod rounds;
#[path = "icprog/rule.rs"]
mod rule;
#[path = "icprog/runs.rs"]
mod runs;
#[path = "icprog/sat_control.rs"]
mod sat_control;
#[path = "icprog/sat_control_publication.rs"]
mod sat_control_publication;
#[path = "icprog/sat_query_law.rs"]
mod sat_query_law;
#[path = "icprog/sat_source.rs"]
mod sat_source;
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

#[derive(Clone, Copy, PartialEq, ValueEnum)]
enum Round {
    R02b,
    R03,
    R05,
    R07,
}

#[derive(Clone, Copy, PartialEq, ValueEnum)]
enum Comparison {
    S23,
    V2,
}

#[derive(Subcommand)]
enum Command {
    /// Replay full published F5 capsule custody as data, without execution.
    F5ControlReplayCustody {
        #[arg(long)]
        publication: PathBuf,
        #[arg(long)]
        registration_sha256: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Publish exact sidecars for an unconsumed F5 capsule and verify its archive.
    F5ControlPublishCustody {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        publication: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Portable data replay of a validation-only F5 build; no archive execution/custody.
    F5ControlReplayValidation {
        #[arg(long)]
        publication: PathBuf,
        #[arg(long)]
        validation_registration_sha256: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Publish compact receipts of a non-executable F5 build-validation capsule.
    F5ControlPublishValidation {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Freeze a disclosed n17 F5 controller/worker/auditor; builds only, no job.
    F5ControlFreeze {
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        cargo: PathBuf,
        #[arg(long)]
        rustc: PathBuf,
        #[arg(long)]
        host_context: PathBuf,
        #[arg(long)]
        validation_only: bool,
    },
    /// Consume the sole invocation of a separately published disclosed F5 seal.
    F5ControlExecute {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
    },
    /// Audit F5 source, transport and arithmetic using its frozen checker; no search.
    F5ControlAudit {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Publish the exact consumed native F5 control; never reruns its worker or audit.
    F5ControlPublishResult {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        audit: PathBuf,
        #[arg(long)]
        source_acceptance: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Portable postexecution F5 data/math replay; executes no archived binary.
    F5ControlReplayResult {
        #[arg(long)]
        publication: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Replay retained disclosed F5 target mathematics; no worker, fresh yield or admission.
    F5TargetReplay {
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Reconstruct retained F5 ordinary preparation natively; no search or new timing.
    F5PreparationReplay {
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Publish the exact consumed native SAT control; creates a new immutable tree.
    SatControlPublish {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        audit: PathBuf,
        #[arg(long)]
        host_context: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Portable postexecution data/math replay; executes no archived binary or solver.
    SatControlReplayPublication {
        #[arg(long)]
        publication: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Freeze disclosed n17 native-control sources, dependencies and binaries; no search.
    SatControlFreeze {
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        config: PathBuf,
        #[arg(long)]
        cargo: PathBuf,
        #[arg(long)]
        rustc: PathBuf,
    },
    /// Consume one frozen disclosed n17 registration before its sole native execution.
    SatControlExecute {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
    },
    /// Independently audit a consumed control; must use its frozen checker binary.
    SatControlAudit {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Independently replay the retained disclosed n17 SAT source witness; no solver or new target.
    SatSourceReplay {
        #[arg(long, default_value = ".")]
        root: PathBuf,
        /// Create an immutable report; refuses to overwrite any existing file.
        #[arg(long)]
        out: Option<PathBuf>,
    },
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
    /// A round's README tables, rendered from its analysis only (Markdown).
    Table {
        round: Round,
        /// The round's `analysis.json`.
        analysis: PathBuf,
    },
    /// An accepted round's baseline entry, appended to `baselines.json`.
    Baseline {
        /// The round's `analysis.json`.
        analysis: PathBuf,
        /// The run tree whose `host.json` (and `host-resumed.json`) it names.
        #[arg(long)]
        runs: PathBuf,
        /// The entry's identity: baseline, round, class, commit, built_from,
        /// base, src_tree, binary_sha256, source, what_changed, measured,
        /// and `previous`, the baseline whose sizes carry the EC1 identities.
        #[arg(long)]
        meta: PathBuf,
        /// The ledger to append to.
        #[arg(long)]
        ledger: PathBuf,
    },
    /// A callgrind profile split by measurement phase, as JSON (the retired
    /// `harness/callgrind_phases.py`).
    CallgrindPhases {
        /// The directory holding the profile and its parts.
        dir: PathBuf,
        /// The `--callgrind-out-file` basename, without a part suffix.
        stem: String,
        /// How many functions to list, by exclusive instructions.
        #[arg(long, default_value_t = 30)]
        top: usize,
    },
    /// R02b's callgrind control from a run tree's `callgrind/`, as JSON.
    CallgrindControl {
        /// The directory holding the arms' `*.phases.json` and `*.workflow.json`.
        dir: PathBuf,
    },
    /// The single-target rule's comparison: `manifest`, `pin`, `size` and
    /// `all` run its rows; `claims` writes every row's manifests, replay and
    /// checked claim; `analyse` prints the figures.
    Rule {
        /// manifest, pin, size, all, reference, claims or analyse.
        step: String,
        /// Which comparison: ledger §23's (frozen) or the one at baseline v2.
        #[arg(long, value_enum, default_value = "s23")]
        comparison: Comparison,
        /// The repository checkout (default: the current directory).
        #[arg(long, default_value = ".")]
        root: PathBuf,
        /// The comparison's directory (default: the comparison's own).
        #[arg(long)]
        dir: Option<PathBuf>,
        /// The run tree (default: the directory's `runs/`).
        #[arg(long)]
        runs: Option<PathBuf>,
        /// Where `claims/` and `manifests/` go (default: the directory).
        #[arg(long)]
        out: Option<PathBuf>,
        /// The repository whose objects hold the binary's build commit
        /// (default: the checkout).
        #[arg(long)]
        git: Option<PathBuf>,
        /// The size a `size` step runs, as `a,n`.
        #[arg(long)]
        size: Option<String>,
        /// The `ic` binary the rows price.
        #[arg(long)]
        ic: Option<PathBuf>,
        /// The commit `ic` was built from, for a fresh manifest.
        #[arg(long)]
        ic_commit: Option<String>,
        /// The isolation tool (default: `isolated_bench` beside this binary).
        #[arg(long)]
        isolate: Option<PathBuf>,
        /// The strong rho fixture a `reference` step walks
        /// (`examples/koblitz_rho_fixture.rs`).
        #[arg(long)]
        fixture: Option<PathBuf>,
        /// The commit the fixture was built from, for a fresh reference check.
        #[arg(long)]
        fixture_commit: Option<String>,
    },
    /// A round's fresh holdouts by suite v1's own construction
    /// (`suite/v1/make_suite.py`, ported): written once into the round's
    /// `holdouts/` with their `SHA256SUMS`, or, with `--check`, re-derived
    /// and compared byte for byte with the files already there.
    Holdouts {
        /// The round's directory, e.g. `research/ic_tool_program/rounds/R07-main-head`.
        round_dir: PathBuf,
        /// The repository checkout (default: the current directory).
        #[arg(long, default_value = ".")]
        root: PathBuf,
        /// A size, as `a,n`; repeat for each.
        #[arg(long = "size", required = true)]
        sizes: Vec<String>,
        /// The recipe seeds, two targets to a seed, in order.
        #[arg(long, value_delimiter = ',', required = true)]
        seeds: Vec<i128>,
        /// The first target's number; the rest follow on.
        #[arg(long)]
        first_target: i128,
        /// Compare with the files already there instead of writing.
        #[arg(long)]
        check: bool,
    },
    /// One of a round's declared steps, natively (R05, R02b, R07).
    Run {
        round: Round,
        /// R05: plan, manifest-resumed, pin, compare, holdout or extend.
        /// R02b and R07: plan, manifest, pin, aa, compare, holdout, extend,
        /// callgrind or manifest-resumed.
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
        /// The commit the base arm was built from, for a fresh manifest.
        #[arg(long)]
        base_commit: Option<String>,
        /// The commit the candidate was built from, for a fresh manifest.
        #[arg(long)]
        cand_commit: Option<String>,
    },
}

/// The programme's directory and the round's.
fn round_dir(round: Round, root: &std::path::Path) -> Result<(PathBuf, PathBuf), String> {
    let programme = suite::programme(root)?;
    let dir = match round {
        Round::R02b => "R02b-wide-tail-retest",
        Round::R03 => "R03-curve-construction",
        Round::R05 => "R05-presence-filter",
        Round::R07 => "R07-main-head",
    };
    let round_dir = programme.join("rounds").join(dir);
    Ok((programme, round_dir))
}

fn context(
    round: Round,
    root: &std::path::Path,
    runs: Option<PathBuf>,
) -> Result<rounds::Ctx, String> {
    let (programme, round_dir) = round_dir(round, root)?;
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

struct RunArgs {
    root: PathBuf,
    runs: Option<PathBuf>,
    arms: [PathBuf; 2],
    isolate: Option<PathBuf>,
    commits: [Option<String>; 2],
}

fn run(round: Round, step: &str, args: RunArgs) -> Result<String, String> {
    let RunArgs {
        root,
        runs,
        arms,
        isolate,
        commits,
    } = args;
    // A fresh round's first step makes its run tree; R03's is frozen.
    let runs = match runs {
        Some(r) => r,
        None => round_dir(round, &root)?.1.join("runs"),
    };
    if round != Round::R03 {
        std::fs::create_dir_all(&runs).map_err(|e| format!("{}: {e}", runs.display()))?;
    }
    let ctx = context(round, &root, Some(runs))?;
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
        Round::R02b => rounds::r02b::run(&ctx, step, &b, &arms, &root, &commits)?,
        Round::R05 => rounds::r05::run(&ctx, step, &b, &arms, &root)?,
        Round::R07 => rounds::r07::run(&ctx, step, &b, &arms, &root, &commits)?,
        Round::R03 => return Err("R03 is complete; its runs are frozen".into()),
    };
    Ok(match doc {
        json::J::Null => String::new(),
        d => json::dumps(&d, 1),
    })
}

/// Append an accepted round's entry to the ledger; none is ever edited.
fn baseline(
    analysis: &std::path::Path,
    runs: &std::path::Path,
    meta: &std::path::Path,
    ledger: &std::path::Path,
) -> Result<String, String> {
    let meta = json::read(meta)?;
    let mut doc = json::read(ledger)?;
    let name = meta.at("baseline")?.clone();
    let host_resumed = json::read_opt(&runs.join("host-resumed.json"))?;
    let entry = report::baseline_entry(
        &json::read(analysis)?,
        &doc,
        &json::read(&runs.join("host.json"))?,
        host_resumed.as_ref(),
        &meta,
    )?;
    let Some(json::J::Arr(list)) = (match &mut doc {
        json::J::Obj(kv) => kv
            .iter_mut()
            .find(|(k, _)| k == "baselines")
            .map(|(_, v)| v),
        _ => None,
    }) else {
        return Err("the ledger has no `baselines` list".into());
    };
    if list.iter().any(|b| b.get("baseline") == Some(&name)) {
        return Err(format!("the ledger already has {}", json::dumps(&name, 0)));
    }
    list.push(entry);
    std::fs::write(ledger, json::dumps_utf8(&doc, 1) + "\n")
        .map_err(|e| format!("{}: {e}", ledger.display()))?;
    Ok(String::new())
}

struct RuleArgs {
    comparison: Comparison,
    root: PathBuf,
    dir: Option<PathBuf>,
    runs: Option<PathBuf>,
    out: Option<PathBuf>,
    git: Option<PathBuf>,
    size: Option<String>,
    ic: Option<PathBuf>,
    ic_commit: Option<String>,
    isolate: Option<PathBuf>,
    fixture: Option<PathBuf>,
    fixture_commit: Option<String>,
}

/// One step of the rule's comparison.  Paths are made absolute but not
/// resolved, so a record's path is relative to the checkout as given.
fn rule(step: &str, args: RuleArgs) -> Result<String, String> {
    let abs =
        |p: &std::path::Path| std::path::absolute(p).map_err(|e| format!("{}: {e}", p.display()));
    let root = abs(&args.root)?;
    let (default_dir, constants) = match args.comparison {
        Comparison::S23 => ("research/ic_single_target_20260930", &rule::S23),
        Comparison::V2 => ("research/ic_tool_program/rule/v2", &rule::V2),
    };
    let here = abs(&args.dir.clone().unwrap_or_else(|| root.join(default_dir)))?;
    let comparison = rule::Comparison {
        git: abs(args.git.as_deref().unwrap_or(&root))?,
        runs: abs(&args.runs.clone().unwrap_or_else(|| here.join("runs")))?,
        out: abs(args.out.as_deref().unwrap_or(&here))?,
        root,
        here,
        constants,
    };
    let tools = || -> Result<rule::Tools, String> {
        let ic = args.ic.clone().ok_or("this step needs --ic")?;
        if !ic.exists() {
            return Err(format!("no binary at {}", ic.display()));
        }
        let isolate = match &args.isolate {
            Some(p) => p.clone(),
            None => std::env::current_exe()
                .map_err(|e| e.to_string())?
                .with_file_name("isolated_bench"),
        };
        if !isolate.exists() {
            return Err(format!("no isolation tool at {}", isolate.display()));
        }
        Ok(rule::Tools {
            ic: abs(&ic)?,
            ic_commit: args.ic_commit.clone(),
            isolate: abs(&isolate)?,
        })
    };
    match step {
        "manifest" => comparison
            .manifest(&tools()?)
            .map(|doc| json::dumps(&doc, 1)),
        "pin" => comparison.pin(&tools()?).map(|doc| json::dumps(&doc, 1)),
        "size" => {
            let size = args.size.as_deref().ok_or("a size step needs --size a,n")?;
            let (a, n) = size
                .split_once(',')
                .and_then(|(a, n)| Some((a.trim().parse().ok()?, n.trim().parse().ok()?)))
                .ok_or("--size is a,n")?;
            if !rule::SIZES.contains(&(a, n)) {
                return Err(format!("k{a}n{n} is not one of the comparison's sizes"));
            }
            comparison
                .size_runs(a, n, &tools()?)
                .map(|()| String::new())
        }
        "all" => comparison.all(&tools()?).map(|()| String::new()),
        "reference" => {
            let fixture = args
                .fixture
                .as_deref()
                .ok_or("a reference step needs --fixture")?;
            if !fixture.exists() {
                return Err(format!("no fixture at {}", fixture.display()));
            }
            let isolate = match &args.isolate {
                Some(p) => p.clone(),
                None => std::env::current_exe()
                    .map_err(|e| e.to_string())?
                    .with_file_name("isolated_bench"),
            };
            comparison
                .reference(
                    &abs(fixture)?,
                    args.fixture_commit.as_deref(),
                    &abs(&isolate)?,
                )
                .map(|()| String::new())
        }
        "claims" => comparison.claims(),
        "analyse" => comparison.analyse().map(|doc| json::dumps(&doc, 1)),
        other => Err(format!(
            "unknown step {other:?}: manifest, pin, size, all, reference, claims or analyse"
        )),
    }
}

fn analyse(round: Round, root: PathBuf, runs: Option<PathBuf>) -> Result<String, String> {
    let ctx = context(round, &root, runs)?;
    let doc = match round {
        Round::R02b => rounds::r02b::analyse(&ctx)?,
        Round::R03 => rounds::r03::analyse(&ctx)?,
        Round::R05 => rounds::r05::analyse(&ctx)?,
        Round::R07 => rounds::r07::analyse(&ctx)?,
    };
    Ok(json::dumps(&doc, 1))
}

struct HoldoutArgs {
    round_dir: PathBuf,
    root: PathBuf,
    sizes: Vec<String>,
    seeds: Vec<i128>,
    first_target: i128,
    check: bool,
}

/// Write a round's fresh holdouts once, or check the ones it has.
fn holdouts(args: HoldoutArgs) -> Result<String, String> {
    let programme = suite::programme(&args.root)?;
    let mut sizes = Vec::new();
    for s in &args.sizes {
        let size = s
            .split_once(',')
            .and_then(|(a, n)| Some((a.trim().parse().ok()?, n.trim().parse().ok()?)))
            .ok_or_else(|| format!("--size {s:?} is not a,n"))?;
        sizes.push(size);
    }
    let files = suite::holdouts(&programme, &sizes, &args.seeds, args.first_target)?;
    let sums = suite::sums_text(&files);
    let dir = args.round_dir.join("holdouts");
    if args.check {
        let mut mismatches = Vec::new();
        for (path, text) in &files {
            if std::fs::read(dir.join(path)).ok().as_deref() != Some(text.as_bytes()) {
                mismatches.push(json::J::Str(path.clone()));
            }
        }
        if std::fs::read_to_string(dir.join("SHA256SUMS"))
            .ok()
            .as_deref()
            != Some(sums.as_str())
        {
            mismatches.push(json::J::Str("SHA256SUMS".into()));
        }
        let doc = json::obj([
            ("files", json::J::Int(files.len() as i128)),
            ("mismatches", json::J::Arr(mismatches.clone())),
        ]);
        return if mismatches.is_empty() {
            Ok(json::dumps(&doc, 1))
        } else {
            Err(json::dumps(&doc, 1))
        };
    }
    if dir.exists() {
        return Err(format!(
            "{} exists; a round's holdouts are drawn once (use --check)",
            dir.display()
        ));
    }
    for (path, text) in &files {
        let out = dir.join(path);
        if let Some(parent) = out.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        std::fs::write(&out, text).map_err(|e| format!("{}: {e}", out.display()))?;
    }
    let out = dir.join("SHA256SUMS");
    std::fs::write(&out, &sums).map_err(|e| format!("{}: {e}", out.display()))?;
    Ok(format!(
        "{} holdout files written to {}",
        files.len(),
        dir.display()
    ))
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::F5ControlReplayCustody {
            publication,
            registration_sha256,
            out,
        } => f5_control::replay_custody(&publication, &registration_sha256, &out),
        Command::F5ControlPublishCustody {
            capsule,
            publication,
            out,
        } => f5_control::publish_custody(&capsule, &publication, &out),
        Command::F5ControlReplayValidation {
            publication,
            validation_registration_sha256,
            out,
        } => f5_control::replay_validation(&publication, &validation_registration_sha256, &out),
        Command::F5ControlPublishValidation { capsule, out } => {
            f5_control::publish_validation(&capsule, &out)
        }
        Command::F5ControlFreeze {
            root,
            out,
            config,
            cargo,
            rustc,
            host_context,
            validation_only,
        } => f5_control::freeze(
            &root,
            &out,
            &config,
            &cargo,
            &rustc,
            &host_context,
            validation_only,
        ),
        Command::F5ControlExecute {
            capsule,
            execution,
            registration_sha256,
        } => f5_control::execute(&capsule, &execution, &registration_sha256),
        Command::F5ControlAudit {
            capsule,
            execution,
            registration_sha256,
            out,
        } => f5_control::audit(&capsule, &execution, &registration_sha256, &out),
        Command::F5ControlPublishResult {
            capsule,
            execution,
            audit,
            source_acceptance,
            out,
        } => {
            f5_control_publication::publish(&capsule, &execution, &audit, &source_acceptance, &out)
        }
        Command::F5ControlReplayResult { publication, out } => {
            f5_control_publication::replay_to(&publication, &out)
        }
        Command::F5TargetReplay { root, out } => f5_target::replay(&root, &out),
        Command::F5PreparationReplay { root, out } => prepared_f5::run(&root, &out),
        Command::SatControlPublish {
            capsule,
            execution,
            audit,
            host_context,
            out,
        } => sat_control_publication::publish(&capsule, &execution, &audit, &host_context, &out),
        Command::SatControlReplayPublication { publication, out } => {
            sat_control_publication::replay_to(&publication, &out)
        }
        Command::SatControlFreeze {
            root,
            out,
            config,
            cargo,
            rustc,
        } => sat_control::freeze(&root, &out, &config, &cargo, &rustc),
        Command::SatControlExecute {
            capsule,
            execution,
            registration_sha256,
        } => sat_control::execute(&capsule, &execution, &registration_sha256),
        Command::SatControlAudit {
            capsule,
            execution,
            out,
        } => sat_control::audit(&capsule, &execution, &out),
        Command::Analyse { round, root, runs } => analyse(round, root, runs),
        Command::SatSourceReplay { root, out } => sat_source::run(&root, out.as_deref()),
        Command::Baseline {
            analysis,
            runs,
            meta,
            ledger,
        } => baseline(&analysis, &runs, &meta, &ledger),
        Command::CallgrindPhases { dir, stem, top } => {
            callgrind::phases(&dir, &stem, top).map(|doc| json::dumps(&doc, 1))
        }
        Command::CallgrindControl { dir } => {
            callgrind::control(&dir).map(|doc| json::dumps(&doc, 1))
        }
        Command::Rule {
            step,
            comparison,
            root,
            dir,
            runs,
            out,
            git,
            size,
            ic,
            ic_commit,
            isolate,
            fixture,
            fixture_commit,
        } => rule(
            &step,
            RuleArgs {
                comparison,
                root,
                dir,
                runs,
                out,
                git,
                size,
                ic,
                ic_commit,
                isolate,
                fixture,
                fixture_commit,
            },
        ),
        Command::Table { round, analysis } => json::read(&analysis).and_then(|doc| match round {
            Round::R05 | Round::R02b | Round::R07 => report::r05(&doc),
            Round::R03 => Err("R03's tables are in its README, written before icprog".into()),
        }),
        Command::Holdouts {
            round_dir,
            root,
            sizes,
            seeds,
            first_target,
            check,
        } => holdouts(HoldoutArgs {
            round_dir,
            root,
            sizes,
            seeds,
            first_target,
            check,
        }),
        Command::Run {
            round,
            step,
            root,
            runs,
            base,
            cand,
            isolate,
            base_commit,
            cand_commit,
        } => run(
            round,
            &step,
            RunArgs {
                root,
                runs,
                arms: [base, cand],
                isolate,
                commits: [base_commit, cand_commit],
            },
        ),
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
