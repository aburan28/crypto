//! `ecbench` — one harness for every ECDLP method.
//!
//! ```text
//! ecbench methods                      what can be measured
//! ecbench host                         the host capsule and its class id
//! ecbench plan   --spec S              validate a spec, print its workloads and identities
//! ecbench run    --spec S --out D      run a session (pinned, interleaved, sealed records)
//! ecbench verify --dir D [--replay N]  re-derive a session; replay runs exactly
//! ecbench compare --dir D --a X --b Y  paired ratio with a bootstrap interval
//! ecbench table  --dir D... [--reference ARM]   the one-unit table AGENTS.md §2 asks for
//! ecbench fb --curve C --factor-base F a factor base with its points
//! ecbench claim build --dir D --ic X --rho Y --workload W   a vs_rho claim, checked
//! ecbench claim attach --base-report R --independent-receipt A   attach a later replay
//! ecbench claim check --report R       the vs_rho checker
//! ecbench db sql PATHS... | sqlite3 ecbench.db
//! ecbench isolab-job --spec S --binary B    an isolab.job/v1 for an independent runner
//! ecbench bound fit --dir D... --arm A --out B.json     a sealed bound: the method's fitted cost
//! ecbench bound check --record B.json      re-derive a committed bound from its sessions
//! ecbench frontier build --bounds DIR --out F.json --markdown FRONTIER.md [--check]
//! ecbench challenge seal --draft C.json --out C.json    validate and seal a challenge
//! ecbench challenge spec --challenge C --candidate '{"id":"..."}' --epoch N --out spec.json
//! ecbench challenge verdict --challenge C --dir D --epoch N --replay-all --out V.json
//! ```
//!
//! Bounds, frontiers and challenges are `docs/bounds/README.md`.
//!
//! The standard is `docs/ecbench/README.md`.

#[path = "ecbench/pdp3_koblitz.rs"]
mod pdp3_koblitz;

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};

use crypto_lib::cryptanalysis::ecbench::{
    audit, bounds, callgrind, canonical, challenge, claim, compare, db, frontier, host, isolab,
    isolation, methods, record, runner, signals, spec, stats, workload,
};

#[derive(Parser)]
#[command(
    name = "ecbench",
    version,
    about = "One harness for every ECDLP method: rho, BSGS, kangaroo, index calculus"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// List every method with its parameters.
    Methods {
        #[arg(long)]
        json: bool,
    },
    /// Print the host capsule.
    Host {
        #[arg(long)]
        json: bool,
    },
    /// Validate a spec and print what it expands to.
    Plan {
        #[arg(long)]
        spec: PathBuf,
        #[arg(long)]
        json: bool,
    },
    /// Run a session.
    Run {
        #[arg(long)]
        spec: PathBuf,
        /// A new directory; `${VAR}` is expanded from the environment.
        #[arg(long)]
        out: String,
        /// auto (a whole core away from CPU 0), inherit (this process's
        /// placement, e.g. isolab's), none (record only, L0), or a list.
        #[arg(long, default_value = "auto")]
        cpus: String,
        #[arg(long, default_value = isolation::DEFAULT_LOCK)]
        lock: String,
        /// Queue behind another benchmark holding the lock.
        #[arg(long)]
        wait: bool,
        /// Record a busy preflight and continue; every run stays below L2.
        #[arg(long)]
        allow_busy: bool,
        #[arg(long, default_value_t = 2.0)]
        settle: f64,
        #[arg(long, default_value_t = 0.10)]
        max_other_cpu: f64,
        #[arg(long, default_value_t = 5.0)]
        max_psi: f64,
        #[arg(long)]
        quiet: bool,
    },
    /// The measured child: reads a job on stdin, prints its output.
    #[command(hide = true)]
    Exec,
    /// Export one exact measured-child input from a frozen spec and sequence.
    ProfileInput {
        #[arg(long)]
        spec: PathBuf,
        #[arg(long)]
        seq: u64,
    },
    /// Sum the Callgrind Ir parts inside the solve markers.
    CallgrindIr {
        /// The value supplied to --callgrind-out-file.
        #[arg(long)]
        prefix: PathBuf,
    },
    /// Split one profiled solve into reusable, one-target online and tail Ir.
    CallgrindOnlineIr {
        /// The value supplied to --callgrind-out-file.
        #[arg(long)]
        prefix: PathBuf,
    },
    /// Audit a session from its own files.
    Verify {
        #[arg(long)]
        dir: String,
        /// Re-execute this many measured runs and require identical counts.
        #[arg(long, default_value_t = 0)]
        replay: usize,
        /// Re-execute every measured deterministic run.
        #[arg(long)]
        replay_all: bool,
        /// Accept a session marked `interrupted`: every integrity check
        /// still applies, but not "every planned execution recorded".
        #[arg(long)]
        allow_interrupted: bool,
        /// Write the receipt here (it is also printed).
        #[arg(long)]
        out: Option<PathBuf>,
        /// Exit 1 when the audit finds a problem (isolab's verify contract).
        #[arg(long)]
        exit_code: bool,
    },
    /// Compare arm B against arm A.
    Compare {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        a: String,
        #[arg(long)]
        b: String,
        /// B's session, when it is another one (operations only).
        #[arg(long)]
        b_dir: Option<PathBuf>,
        #[arg(long, default_value_t = 4000)]
        resamples: usize,
        #[arg(long, default_value_t = 20261002)]
        seed: u64,
        /// Save to <dir>/comparisons/<a>__<b>.json.
        #[arg(long)]
        save: bool,
        #[arg(long)]
        json: bool,
    },
    /// The one-unit table: every arm on every curve, mean S, ratio to the
    /// floor and to the reference arm.
    Table {
        #[arg(long, num_args = 1..)]
        dir: Vec<PathBuf>,
        /// The arm every other is divided by (default: the session's
        /// `reference` arm).
        #[arg(long)]
        reference: Option<String>,
    },
    /// Build a factor base and write it with its points.
    Fb {
        /// A curve spec, e.g. '{"kind":"koblitz","a":0,"n":31}'.
        #[arg(long)]
        curve: String,
        /// A factor-base plug-in spec, e.g. 'koblitz-orbit:divisor=1'.
        #[arg(long)]
        factor_base: String,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// `vs_rho` claims: one IC run against the strong rho run on one
    /// public target, in the report `docs/ic/boundary_targets.json` defines.
    Claim {
        #[command(subcommand)]
        cmd: ClaimCmd,
    },
    /// The database.
    Db {
        #[command(subcommand)]
        cmd: DbCmd,
    },
    /// An isolab.job/v1 that runs a spec on an independent lab worker.
    IsolabJob {
        #[arg(long)]
        spec: PathBuf,
        /// Upload this binary (built for the worker's architecture).
        #[arg(long, conflicts_with_all = ["git_url", "commit"])]
        binary: Option<String>,
        /// Or build from this repository at `--commit` on the worker.
        #[arg(long, requires = "commit")]
        git_url: Option<String>,
        #[arg(long)]
        commit: Option<String>,
        #[arg(long, default_value = "default")]
        pool: String,
        #[arg(long, default_value_t = 2)]
        cpus: u32,
        #[arg(long)]
        memory_mb: Option<u64>,
        #[arg(long)]
        image: Option<String>,
        #[arg(long, default_value = "strict")]
        policy: String,
        #[arg(long, default_value_t = 0)]
        timeout_seconds: u64,
    },
    /// Bound records: what a method costs on a domain, fitted across
    /// sizes, sealed and re-derivable (`docs/bounds/README.md`).
    Bound {
        #[command(subcommand)]
        cmd: BoundCmd,
    },
    /// The frontier of every domain's bounds, and its page.
    Frontier {
        #[command(subcommand)]
        cmd: FrontierCmd,
    },
    /// Challenges: the frozen paired spec a candidate runs against an
    /// incumbent, and the verdict that moves (or does not move) the frontier.
    Challenge {
        #[command(subcommand)]
        cmd: ChallengeCmd,
    },
}

#[derive(Subcommand)]
enum BoundCmd {
    /// Fit the bound of one arm from one or more sessions of it.
    Fit {
        #[arg(long, num_args = 1..)]
        dir: Vec<PathBuf>,
        #[arg(long)]
        arm: String,
        /// Keep only curves of this tier (toy, medium, crypto).
        #[arg(long)]
        tier: Option<String>,
        /// Keep only these curve slugs.
        #[arg(long)]
        curve: Vec<String>,
        /// Audit receipts of the sessions, cited in the provenance.
        #[arg(long)]
        audit: Vec<PathBuf>,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "")]
        notes: String,
        /// Paths in the record are written relative to this directory.
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long, default_value_t = bounds::RESAMPLES)]
        resamples: usize,
        #[arg(long, default_value_t = bounds::SEED)]
        seed: u64,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Re-derive committed bounds from the sessions they name and require
    /// the same ids.
    Check {
        #[arg(long, num_args = 1..)]
        record: Vec<PathBuf>,
        #[arg(long, default_value = ".")]
        root: PathBuf,
    },
}

#[derive(Subcommand)]
enum FrontierCmd {
    /// Build the frontier from bound records (files or directories).
    Build {
        #[arg(long, num_args = 1..)]
        bounds: Vec<PathBuf>,
        /// Dominance axes, comma-separated: ops, memory, uncharged,
        /// field_muls, field_sqrs, field_invs.
        #[arg(long, default_value = "ops,memory")]
        axes: String,
        #[arg(long)]
        out: Option<PathBuf>,
        #[arg(long)]
        markdown: Option<PathBuf>,
        /// Write nothing; fail if --out or --markdown is stale.
        #[arg(long)]
        check: bool,
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand)]
enum ChallengeCmd {
    /// Validate a drafted challenge (its curves build, lie in its domain,
    /// and number enough sizes) and seal it.
    Seal {
        #[arg(long)]
        draft: PathBuf,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Check committed challenges: seals, domain ids, curves.
    Check {
        #[arg(long, num_args = 1..)]
        file: Vec<PathBuf>,
    },
    /// The spec a candidate runs against the challenge in an epoch.
    Spec {
        #[arg(long)]
        challenge: PathBuf,
        /// The candidate as a method spec, e.g. '{"id":"rho.plain"}'.
        #[arg(long)]
        candidate: String,
        #[arg(long)]
        epoch: u64,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Judge a session run from `challenge spec`: audit, pair, fit, decide.
    Verdict {
        #[arg(long)]
        challenge: PathBuf,
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        epoch: u64,
        #[arg(long, default_value = challenge::CANDIDATE_ARM)]
        candidate: String,
        #[arg(long, default_value = challenge::INCUMBENT_ARM)]
        incumbent: String,
        #[arg(long, default_value = challenge::CONTROL_ARM)]
        control: String,
        /// Runs to replay in the audit.
        #[arg(long, default_value_t = 0)]
        replay: usize,
        /// Replay every measured deterministic run (what acceptance
        /// requires by default).
        #[arg(long)]
        replay_all: bool,
        /// Committed bound records, to compare the incumbent's recorded
        /// figure with what it measured here.
        #[arg(long)]
        bounds: Option<PathBuf>,
        #[arg(long, default_value = ".")]
        root: PathBuf,
        #[arg(long, default_value_t = bounds::RESAMPLES)]
        resamples: usize,
        #[arg(long, default_value_t = bounds::SEED)]
        seed: u64,
        #[arg(long)]
        out: Option<PathBuf>,
        /// Write the candidate's new bound here on an advance or trade.
        #[arg(long)]
        bound_out: Option<PathBuf>,
        /// Write the audit receipt here.
        #[arg(long)]
        audit_out: Option<PathBuf>,
        /// Exit 1 when the verdict is inadmissible.
        #[arg(long)]
        exit_code: bool,
    },
}

#[derive(Subcommand)]
enum ClaimCmd {
    /// Assemble the claim for one workload and round of a session, and
    /// check it.  Run it with the binary that measured the session.
    Build {
        #[arg(long)]
        dir: PathBuf,
        /// The index-calculus arm.
        #[arg(long)]
        ic: String,
        /// The `rho.signed_frobenius_strong` arm.
        #[arg(long)]
        rho: String,
        /// The ecbench workload id (`W` + 12 hex).
        #[arg(long)]
        workload: String,
        /// The round, counted from 0 with warm-ups included (default: the
        /// first measured round both arms ran).
        #[arg(long)]
        round: Option<u32>,
        /// An audit receipt (`ecbench verify --replay-all --out`) made on
        /// another host class that reproduced both runs.
        #[arg(long, requires = "pointer")]
        independent_receipt: Option<PathBuf>,
        /// Where that receipt is durably stored.
        #[arg(long, requires = "independent_receipt")]
        pointer: Option<String>,
        /// Write the claim here (it is printed otherwise).
        #[arg(long)]
        out: Option<PathBuf>,
        /// Exit 1 when the checker finds a problem.
        #[arg(long)]
        exit_code: bool,
    },
    /// Attach a later other-host replay to an unchanged claim produced by
    /// the measuring binary. The entire original report is reconstructed
    /// from the frozen session before provenance is added.
    Attach {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        ic: String,
        #[arg(long)]
        rho: String,
        #[arg(long)]
        workload: String,
        #[arg(long)]
        round: Option<u32>,
        #[arg(long)]
        base_report: PathBuf,
        #[arg(long)]
        independent_receipt: PathBuf,
        #[arg(long)]
        pointer: String,
        #[arg(long)]
        out: Option<PathBuf>,
        #[arg(long)]
        exit_code: bool,
    },
    /// Check a report against the `vs_rho` schema, as
    /// `boundary_autolab.py claim-check --stage vs_rho` does; exit 1 on
    /// any problem.
    Check {
        #[arg(long)]
        report: PathBuf,
    },
}

#[derive(Subcommand)]
enum DbCmd {
    /// Print the schema.
    Schema,
    /// Print SQL for session directories, comparison files and factor-base dumps.
    Sql {
        #[arg(required = true)]
        paths: Vec<PathBuf>,
    },
}

fn main() -> ExitCode {
    // Oracles this binary links that the library cannot name (see
    // `pdp3_koblitz`); registered before any run, child or replay.
    crypto_lib::cryptanalysis::ecbench::methods::register_binary_plugins(pdp3_koblitz::plugins());
    match run(Cli::parse()) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("ecbench: {e}");
            ExitCode::from(2)
        }
    }
}

fn read(p: &Path) -> Result<String, String> {
    std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))
}

fn report_problems(problems: &[String]) {
    if problems.is_empty() {
        eprintln!("vs_rho check: passes");
    } else {
        eprintln!("vs_rho check: {} problem(s)", problems.len());
        for p in problems {
            eprintln!("  - {p}");
        }
    }
}

fn print_json(v: &impl serde::Serialize) -> Result<(), String> {
    println!(
        "{}",
        serde_json::to_string_pretty(v).map_err(|e| e.to_string())?
    );
    Ok(())
}

fn load_independent(path: &Path, pointer: String) -> Result<claim::Independent, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(claim::Independent {
        receipt: serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))?,
        receipt_sha256: canonical::sha256_hex(&bytes),
        pointer,
    })
}

fn emit_claim(
    report: serde_json::Value,
    out: Option<PathBuf>,
    exit_code: bool,
) -> Result<ExitCode, String> {
    let problems = claim::check_problems(&claim::check_vs_rho(&report)?);
    let text = serde_json::to_string_pretty(&report).map_err(|e| e.to_string())? + "\n";
    match out {
        Some(o) => {
            std::fs::write(&o, &text).map_err(|e| format!("{}: {e}", o.display()))?;
            eprintln!(
                "claim {} sha256 {}",
                o.display(),
                canonical::sha256_hex(text.as_bytes())
            );
        }
        None => print!("{text}"),
    }
    eprintln!("{}", report["verdict"].as_str().unwrap_or(""));
    report_problems(&problems);
    if exit_code && !problems.is_empty() {
        return Ok(ExitCode::from(1));
    }
    Ok(ExitCode::SUCCESS)
}

fn run(cli: Cli) -> Result<ExitCode, String> {
    match cli.cmd {
        Cmd::Methods { json } => {
            if json {
                let v: Vec<_> = methods::registry()
                    .iter()
                    .map(|m| {
                        serde_json::json!({
                            "id": m.id, "family": m.family, "summary": m.summary, "entry": m.entry,
                            "applies": m.applies,
                            "params": m.params.iter().map(|p| serde_json::json!({
                                "name": p.name, "default": p.default, "help": p.help})).collect::<Vec<_>>(),
                        })
                    })
                    .collect();
                print_json(&v)?;
            } else {
                for m in methods::registry() {
                    println!("{:<22} {:<9} {}", m.id, m.family, m.summary);
                    for p in m.params {
                        println!(
                            "    {:<22} {:<18} {}",
                            p.name,
                            p.default
                                .map(|d| format!("default {d:?}"))
                                .unwrap_or_else(|| "REQUIRED".into()),
                            p.help
                        );
                    }
                }
            }
        }
        Cmd::Host { json } => {
            let c = host::capture()?;
            if json {
                print_json(&c)?;
            } else {
                let s = &c.stable;
                println!("env class   {}", c.env_class_id);
                println!(
                    "cpu         {} ({} logical, {:?} cores, {} NUMA nodes)",
                    s.cpu_model.as_deref().unwrap_or("?"),
                    s.logical_cpus,
                    s.physical_cores,
                    s.numa_nodes
                );
                println!("features    {}", s.features.join(" "));
                println!(
                    "os          {} {} {}",
                    s.os,
                    s.arch,
                    s.kernel_release.as_deref().unwrap_or("")
                );
                println!("virtual     {:?}", host::is_virtual(s));
                println!(
                    "binary      {}",
                    c.build.binary_sha256.as_deref().unwrap_or("?")
                );
                println!(
                    "commit      {} (dirty: {:?}, source: {})",
                    c.build.git_commit.as_deref().unwrap_or("?"),
                    c.build.git_dirty,
                    c.build.git_commit_source.as_deref().unwrap_or("?")
                );
                if c.build.runtime_git_commit != c.build.git_commit {
                    println!(
                        "runtime git {} (dirty: {:?}; diagnostic only)",
                        c.build.runtime_git_commit.as_deref().unwrap_or("none"),
                        c.build.runtime_git_dirty
                    );
                }
                if c.build.debug_assertions {
                    println!(
                        "note        debug build: operation counts are valid, wall time is not"
                    );
                }
                if !s.perf_levels.is_empty() {
                    println!(
                        "note        heterogeneous cores {:?}: no affinity control, runs earn L0",
                        s.perf_levels.iter().map(|l| &l.name).collect::<Vec<_>>()
                    );
                }
            }
        }
        Cmd::Plan { spec: path, json } => {
            let text = read(&path)?;
            let p = spec::plan(spec::Spec::from_json(&text)?)?;
            if json {
                print_json(&serde_json::json!({
                    "spec_id": p.spec_id, "spec_sha256": p.spec_sha256,
                    "arms": p.arms, "workloads": p.workloads, "executions": p.executions.len(),
                }))?;
            } else {
                println!(
                    "spec {} ({}): {} executions",
                    p.spec_id,
                    p.spec.label,
                    p.executions.len()
                );
                for a in &p.arms {
                    println!(
                        "  arm {:<16} {:<10} {} {} {:?}",
                        a.name,
                        format!("{:?}", a.role).to_lowercase(),
                        a.method.method_id,
                        a.method.id,
                        a.method.params
                    );
                }
                for w in &p.workloads {
                    println!(
                        "  workload {} {} r=2^{:.2} A={} registered={} ec1={}",
                        w.workload_id,
                        w.curve.slug,
                        (w.curve.r as f64).log2(),
                        w.curve.automorphisms_available,
                        w.curve.registered,
                        w.curve.ec1.as_deref().unwrap_or("-")
                    );
                }
                let unregistered: Vec<_> = p
                    .workloads
                    .iter()
                    .filter(|w| !w.curve.registered)
                    .map(|w| &w.curve.slug)
                    .collect();
                if !unregistered.is_empty() {
                    println!("note: {} slug(s) are not in docs/curves/registry.json; register before citing them in prose (AGENTS.md §11)", unregistered.len());
                }
            }
        }
        Cmd::Run {
            spec: path,
            out,
            cpus,
            lock,
            wait,
            allow_busy,
            settle,
            max_other_cpu,
            max_psi,
            quiet,
        } => {
            let text = read(&path)?;
            let out = PathBuf::from(isolab::expand_env(&out)?);
            let mut req = isolation::CpuRequest::parse(&cpus)?;
            if req == isolation::CpuRequest::Auto && isolation::topology().cpus.is_empty() {
                eprintln!("ecbench: no CPU topology on this OS; running with --cpus none (every run is recorded at L0, operation counts unaffected)");
                req = isolation::CpuRequest::None;
            }
            let opts = runner::RunOptions {
                cpus: req,
                lock_path: lock,
                wait_for_lock: wait,
                allow_busy,
                thresholds: isolation::Thresholds {
                    settle_seconds: settle,
                    max_other_cpu,
                    max_psi,
                    ..Default::default()
                },
                // On Linux each child execs this very binary through
                // /proc/self/exe, so a rebuild mid-session cannot swap it.
                exe: if cfg!(target_os = "linux") {
                    PathBuf::from("/proc/self/exe")
                } else {
                    std::env::current_exe().map_err(|e| e.to_string())?
                },
            };
            let result = runner::run_session(&text, &out, opts, |line| {
                if !quiet {
                    // A closed stderr must not panic a running session.
                    let _ = writeln!(std::io::stderr(), "{line}");
                }
            });
            let s = match result {
                Ok(s) => s,
                Err(e) => {
                    let _ = writeln!(std::io::stderr(), "ecbench: {e}");
                    // 128 + signal, the shell's convention, when interrupted.
                    let code = signals::interrupted().map(|sig| 128 + sig).unwrap_or(2);
                    return Ok(ExitCode::from(code.clamp(1, 255) as u8));
                }
            };
            let _ = writeln!(
                std::io::stderr(),
                "session {} {} -> {} ({} records: {:?})",
                s.session_id,
                s.status,
                out.display(),
                s.records_written,
                s.status_counts
            );
        }
        Cmd::Exec => {
            let mut input = String::new();
            std::io::Read::read_to_string(&mut std::io::stdin(), &mut input)
                .map_err(|e| e.to_string())?;
            let out = record::child_main(&input);
            println!(
                "{}",
                serde_json::to_string(&out).map_err(|e| e.to_string())?
            );
        }
        Cmd::ProfileInput { spec: path, seq } => {
            let source =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let p = spec::plan(spec::Spec::from_json(&source)?)?;
            let ex = p
                .executions
                .get(usize::try_from(seq).map_err(|_| "sequence does not fit usize")?)
                .filter(|ex| ex.seq == seq)
                .ok_or_else(|| format!("sequence {seq} is not in the plan"))?;
            let arm = &p.arms[ex.arm];
            let w = &p.workloads[ex.workload];
            let input = record::ChildInput {
                schema: record::CHILD_INPUT_SCHEMA.into(),
                curve: workload::CurveSpec::explicit(p.instance(ex.workload))
                    .unwrap_or_else(|| w.curve_spec.clone()),
                target_seed: w.target_seed,
                target_index: w.target_index,
                target_kind: w.kind(),
                expected_workload_id: w.workload_id.clone(),
                method: methods::MethodSpec {
                    id: arm.method.id.clone(),
                    params: arm.method.params.clone(),
                },
                expected_method_id: arm.method.method_id.clone(),
                algorithm_seed: ex.algorithm_seed,
            };
            print_json(&input)?;
        }
        Cmd::CallgrindIr { prefix } => {
            print_json(&callgrind::solve_ir(&prefix)?)?;
        }
        Cmd::CallgrindOnlineIr { prefix } => {
            print_json(&callgrind::online_ir(&prefix)?)?;
        }
        Cmd::Verify {
            dir,
            replay,
            replay_all,
            allow_interrupted,
            out,
            exit_code,
        } => {
            let dir = PathBuf::from(isolab::expand_env(&dir)?);
            let rep = audit::audit_with(
                &dir,
                if replay_all {
                    audit::REPLAY_ALL
                } else {
                    replay
                },
                allow_interrupted,
            )?;
            let text = serde_json::to_string_pretty(&rep).map_err(|e| e.to_string())? + "\n";
            if let Some(o) = out {
                std::fs::write(&o, &text).map_err(|e| format!("{}: {e}", o.display()))?;
                eprintln!("receipt sha256 {}", canonical::sha256_hex(text.as_bytes()));
            }
            print!("{text}");
            eprintln!(
                "audit {}: {} records, {} verified, {} replays reproduced, {} problems",
                if rep.ok { "OK" } else { "FAILED" },
                rep.records,
                rep.verified_records,
                rep.replays.iter().filter(|r| r.reproduced).count(),
                rep.problems.len()
            );
            if exit_code && !rep.ok {
                return Ok(ExitCode::from(1));
            }
        }
        Cmd::Compare {
            dir,
            a,
            b,
            b_dir,
            resamples,
            seed,
            save,
            json,
        } => {
            let c = compare::compare(&dir, &a, b_dir.as_deref(), &b, resamples, seed)?;
            if save {
                let d = dir.join("comparisons");
                std::fs::create_dir_all(&d).map_err(|e| e.to_string())?;
                let f = d.join(format!("{a}__{b}.json"));
                std::fs::write(
                    &f,
                    serde_json::to_string_pretty(&c).map_err(|e| e.to_string())? + "\n",
                )
                .map_err(|e| e.to_string())?;
                eprintln!("saved {}", f.display());
            }
            if json {
                print_json(&c)?;
            } else {
                println!("{}", c.verdict);
                let o = |v: Option<f64>| v.map(|x| format!("{x:.4}")).unwrap_or_else(|| "–".into());
                let ci = |v: Option<(f64, f64)>| {
                    v.map(|(l, h)| format!("[{l:.4}, {h:.4}]"))
                        .unwrap_or_else(|| "–".into())
                };
                println!("| curve | log2 r | workloads | {b} / {a} | 95% interval | interval |");
                println!("|---|---:|---:|---:|---|---|");
                for k in &c.curves {
                    println!(
                        "| {} | {:.2} | {} | {} | {} | {} |",
                        k.slug,
                        k.log2_r,
                        k.workloads,
                        o(k.ratio_b_over_a),
                        ci(k.ci95),
                        k.ci_method
                    );
                }
                println!();
                println!("| workload | curve | log2 r | {a} mean S | {b} mean S | ratio |");
                println!("|---|---|---:|---:|---:|---:|");
                for w in &c.workloads {
                    println!(
                        "| {} | {} | {:.2} | {} | {} | {} |",
                        w.workload_id,
                        w.slug,
                        (w.r as f64).log2(),
                        o(w.a_mean_s),
                        o(w.b_mean_s),
                        o(w.ratio_b_over_a)
                    );
                }
            }
        }
        Cmd::Table { dir, reference } => table(&dir, reference.as_deref())?,
        Cmd::Fb {
            curve,
            factor_base,
            out,
        } => {
            let spec: workload::CurveSpec =
                serde_json::from_str(&curve).map_err(|e| format!("--curve: {e}"))?;
            let dump = methods::dump_factor_base(&spec, &factor_base)?;
            let text = serde_json::to_string_pretty(&dump).map_err(|e| e.to_string())? + "\n";
            match out {
                Some(p) => {
                    std::fs::write(&p, &text).map_err(|e| e.to_string())?;
                    eprintln!(
                        "{} on {}: {} points, {} columns -> {}",
                        dump.factor_base.fb_id,
                        dump.curve.slug,
                        dump.factor_base.signed_points,
                        dump.factor_base.columns,
                        p.display()
                    );
                }
                None => print!("{text}"),
            }
        }
        Cmd::Claim { cmd } => match cmd {
            ClaimCmd::Build {
                dir,
                ic,
                rho,
                workload,
                round,
                independent_receipt,
                pointer,
                out,
                exit_code,
            } => {
                let independent = match (independent_receipt, pointer) {
                    (Some(p), Some(pointer)) => Some(load_independent(&p, pointer)?),
                    _ => None,
                };
                let report = claim::build(&dir, &ic, &rho, &workload, round, independent.as_ref())?;
                return emit_claim(report, out, exit_code);
            }
            ClaimCmd::Attach {
                dir,
                ic,
                rho,
                workload,
                round,
                base_report,
                independent_receipt,
                pointer,
                out,
                exit_code,
            } => {
                let original: serde_json::Value = serde_json::from_str(&read(&base_report)?)
                    .map_err(|e| format!("{}: {e}", base_report.display()))?;
                let independent = load_independent(&independent_receipt, pointer)?;
                let report = claim::attach_independent(
                    &dir,
                    &ic,
                    &rho,
                    &workload,
                    round,
                    &original,
                    &independent,
                )?;
                return emit_claim(report, out, exit_code);
            }
            ClaimCmd::Check { report } => {
                let v: serde_json::Value = serde_json::from_str(&read(&report)?)
                    .map_err(|e| format!("{}: {e}", report.display()))?;
                let check = claim::check_vs_rho(&v)?;
                print_json(&check)?;
                let problems = claim::check_problems(&check);
                report_problems(&problems);
                if !problems.is_empty() {
                    return Ok(ExitCode::from(1));
                }
            }
        },
        Cmd::Db { cmd } => match cmd {
            DbCmd::Schema => print!("{}", db::SCHEMA_SQL),
            DbCmd::Sql { paths } => {
                let refs: Vec<&Path> = paths.iter().map(|p| p.as_path()).collect();
                print!("{}", db::sql(&refs)?);
            }
        },
        Cmd::IsolabJob {
            spec: path,
            binary,
            git_url,
            commit,
            pool,
            cpus,
            memory_mb,
            image,
            policy,
            timeout_seconds,
        } => {
            let text = read(&path)?;
            let p = spec::plan(spec::Spec::from_json(&text)?)?;
            let delivery = match (binary, git_url, commit) {
                (Some(b), None, None) => isolab::Delivery::Binary { path: b },
                (None, Some(url), Some(commit)) => isolab::Delivery::Source { url, commit },
                _ => return Err("pass --binary PATH, or --git-url URL --commit SHA".into()),
            };
            let job = isolab::job(
                &p,
                &text,
                &isolab::JobOptions {
                    pool,
                    delivery,
                    cpus,
                    memory_mb,
                    image,
                    timeout_seconds,
                    policy,
                },
            )?;
            print_json(&job)?;
        }
        Cmd::Bound { cmd } => match cmd {
            BoundCmd::Fit {
                dir,
                arm,
                tier,
                curve,
                audit,
                label,
                notes,
                root,
                resamples,
                seed,
                out,
            } => {
                let b = bounds::bound_from_sessions(&bounds::FitInputs {
                    dirs: &dir,
                    arm: &arm,
                    audits: &audit,
                    label: label.as_deref(),
                    notes: &notes,
                    root: &root,
                    options: bounds::FitOptions {
                        tier,
                        curves: curve,
                        resamples,
                        seed,
                    },
                })?;
                let mut b = b;
                let text = b.seal()?;
                write_or_print(out.as_deref(), &text)?;
                eprintln!(
                    "{} {} ({}): {} on {} size(s), {} verified runs; ops {:.3}× floor, α {} ; {}",
                    b.bound_id,
                    b.method.id,
                    b.domain.tier,
                    b.level,
                    b.fit.sizes,
                    b.provenance.verified,
                    b.constant.ratio_to_floor.value,
                    b.fit
                        .alpha
                        .map(|a| format!("{a:.3}"))
                        .unwrap_or_else(|| "–".into()),
                    b.admissibility.status
                );
                for r in &b.admissibility.reasons {
                    eprintln!("  - {r}");
                }
            }
            BoundCmd::Check { record, root } => {
                let mut failed = 0;
                for r in &record {
                    match bounds::check(r, &root) {
                        Ok(b) => eprintln!("{}: re-derived {} OK", r.display(), b.bound_id),
                        Err(e) => {
                            failed += 1;
                            eprintln!("{e}");
                        }
                    }
                }
                if failed > 0 {
                    return Ok(ExitCode::from(1));
                }
            }
        },
        Cmd::Frontier { cmd } => match cmd {
            FrontierCmd::Build {
                bounds: paths,
                axes,
                out,
                markdown,
                check,
                json,
            } => {
                let loaded = frontier::load_bounds(&paths)?;
                let axes: Vec<String> = axes
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
                let mut f = frontier::build(&loaded, &axes)?;
                let text = f.seal()?;
                let records_dir = paths
                    .first()
                    .map(|p| p.display().to_string())
                    .unwrap_or_default();
                let page = frontier::render_markdown(&f, &records_dir);
                if check {
                    let mut stale = Vec::new();
                    if let Some(o) = &out {
                        if read(o)? != text {
                            stale.push(o.display().to_string());
                        }
                    }
                    if let Some(m) = &markdown {
                        if read(m)? != page {
                            stale.push(m.display().to_string());
                        }
                    }
                    if !stale.is_empty() {
                        eprintln!(
                            "stale: {}; rebuild with `ecbench frontier build`",
                            stale.join(", ")
                        );
                        return Ok(ExitCode::from(1));
                    }
                    eprintln!("frontier {} is current", f.frontier_id);
                } else {
                    if let Some(o) = &out {
                        std::fs::write(o, &text).map_err(|e| format!("{}: {e}", o.display()))?;
                    }
                    if let Some(m) = &markdown {
                        std::fs::write(m, &page).map_err(|e| format!("{}: {e}", m.display()))?;
                    }
                    if json || (out.is_none() && markdown.is_none()) {
                        print!("{text}");
                    }
                }
                for d in &f.domains {
                    eprintln!(
                        "{} {} {} {}: {} bound(s), {} on the frontier; ops leader {}",
                        d.domain.family,
                        d.domain.target_kind,
                        d.domain.tier,
                        d.domain_id,
                        d.entries.len(),
                        d.entries.iter().filter(|e| e.is_frontier).count(),
                        d.ops_leader.as_deref().unwrap_or("–")
                    );
                }
                if !f.inadmissible.is_empty() {
                    eprintln!("{} inadmissible bound(s) listed", f.inadmissible.len());
                }
            }
        },
        Cmd::Challenge { cmd } => match cmd {
            ChallengeCmd::Seal { draft, out } => {
                let mut v: serde_json::Value = serde_json::from_str(&read(&draft)?)
                    .map_err(|e| format!("{}: {e}", draft.display()))?;
                let domain: bounds::Domain = serde_json::from_value(
                    v.get("domain").cloned().ok_or("the draft has no domain")?,
                )
                .map_err(|e| format!("domain: {e}"))?;
                v["domain_id"] = serde_json::Value::String(domain.id()?);
                v["challenge_id"] = serde_json::Value::String(String::new());
                v["schema"] = serde_json::Value::String(challenge::CHALLENGE_SCHEMA.into());
                let mut c: challenge::Challenge =
                    serde_json::from_value(v).map_err(|e| format!("{}: {e}", draft.display()))?;
                challenge::validate(&c)?;
                let text = c.seal()?;
                write_or_print(out.as_deref(), &text)?;
                eprintln!(
                    "challenge {} ({}): {} on {} curve(s), incumbent {}",
                    c.challenge_id,
                    c.domain.tier,
                    c.domain.family,
                    c.workloads.curves.len(),
                    c.incumbent.method.id
                );
            }
            ChallengeCmd::Check { file } => {
                let mut failed = 0;
                for f in &file {
                    match challenge::Challenge::read(f)
                        .and_then(|c| challenge::validate(&c).map(|_| c))
                    {
                        Ok(c) => eprintln!(
                            "{}: challenge {} OK ({} {} {}, {} curves, incumbent {})",
                            f.display(),
                            c.challenge_id,
                            c.domain.family,
                            c.domain.target_kind,
                            c.domain.tier,
                            c.workloads.curves.len(),
                            c.incumbent.method.id
                        ),
                        Err(e) => {
                            failed += 1;
                            eprintln!("{e}");
                        }
                    }
                }
                if failed > 0 {
                    return Ok(ExitCode::from(1));
                }
            }
            ChallengeCmd::Spec {
                challenge: path,
                candidate,
                epoch,
                out,
            } => {
                let c = challenge::Challenge::read(&path)?;
                challenge::validate(&c)?;
                let m: methods::MethodSpec =
                    serde_json::from_str(&candidate).map_err(|e| format!("--candidate: {e}"))?;
                let s = challenge::spec_for(&c, &m, epoch);
                let p = spec::plan(s.clone())?;
                let text = serde_json::to_string_pretty(&s).map_err(|e| e.to_string())? + "\n";
                write_or_print(out.as_deref(), &text)?;
                eprintln!(
                    "spec {} for challenge {} epoch {epoch}: {} executions, target seed {}",
                    p.spec_id,
                    c.challenge_id,
                    p.executions.len(),
                    s.workloads.target_seed
                );
            }
            ChallengeCmd::Verdict {
                challenge: path,
                dir,
                epoch,
                candidate,
                incumbent,
                control,
                replay,
                replay_all,
                bounds: bdir,
                root,
                resamples,
                seed,
                out,
                bound_out,
                audit_out,
                exit_code,
            } => {
                let c = challenge::Challenge::read(&path)?;
                let o = challenge::verdict(&challenge::VerdictInputs {
                    challenge: &c,
                    dir: &dir,
                    epoch,
                    candidate_arm: &candidate,
                    incumbent_arm: &incumbent,
                    control_arm: Some(&control),
                    replay: if replay_all {
                        audit::REPLAY_ALL
                    } else {
                        replay
                    },
                    bounds: bdir.as_deref(),
                    root: &root,
                    resamples,
                    seed,
                })?;
                if let Some(a) = &audit_out {
                    std::fs::write(a, &o.audit_text)
                        .map_err(|e| format!("{}: {e}", a.display()))?;
                }
                write_or_print(out.as_deref(), &o.verdict_text)?;
                if let (Some(b), Some(text)) = (&bound_out, &o.candidate_bound_text) {
                    std::fs::write(b, text).map_err(|e| format!("{}: {e}", b.display()))?;
                    eprintln!(
                        "new bound {} -> {}",
                        o.candidate_bound
                            .as_ref()
                            .map(|x| x.bound_id.as_str())
                            .unwrap_or(""),
                        b.display()
                    );
                } else if bound_out.is_some() {
                    eprintln!("no new bound: the verdict is {}", o.verdict.outcome);
                }
                eprintln!("verdict {}: {}", o.verdict.verdict_id, o.verdict.statement);
                if exit_code && o.verdict.outcome == "inadmissible" {
                    return Ok(ExitCode::from(1));
                }
            }
        },
    }
    Ok(ExitCode::SUCCESS)
}

fn write_or_print(out: Option<&Path>, text: &str) -> Result<(), String> {
    match out {
        Some(o) => {
            std::fs::write(o, text).map_err(|e| format!("{}: {e}", o.display()))?;
            eprintln!(
                "wrote {} (sha256 {})",
                o.display(),
                canonical::sha256_hex(text.as_bytes())
            );
        }
        None => print!("{text}"),
    }
    Ok(())
}

/// One row per (session, arm, curve): mean S over verified measured runs
/// with its two-stage bootstrap interval (workloads, then runs), the
/// method's derived expectation, and the ratios to the curve's generic
/// floor and to the reference arm on the same workloads.
fn table(dirs: &[PathBuf], reference: Option<&str>) -> Result<(), String> {
    println!("| session | arm | method | curve | log2 r | verified | mean S | 95% interval | theory S | S / theory | S / floor | S / reference | lower bound | levels |");
    println!("|---|---|---|---|---:|---:|---:|---|---:|---:|---:|---:|---|---|");
    for d in dirs {
        let s = runner::read_session(d)?;
        let recs = runner::read_records(d)?;
        let plan: runner::PlanDoc =
            serde_json::from_str(&read(&d.join("plan.json"))?).map_err(|e| e.to_string())?;
        let ref_arm = reference.map(String::from).or_else(|| {
            plan.arms
                .iter()
                .find(|a| a.role == spec::Role::Reference)
                .map(|a| a.name.clone())
        });
        let mut curves: Vec<&str> = plan
            .workloads
            .iter()
            .map(|w| w.curve.slug.as_str())
            .collect();
        curves.dedup();
        for arm in &plan.arms {
            for slug in &curves {
                let mine: Vec<&record::Record> = recs
                    .iter()
                    .filter(|r| !r.warmup && r.arm == arm.name && r.workload.curve.slug == *slug)
                    .collect();
                let ok: Vec<f64> = mine
                    .iter()
                    .filter(|r| r.counts())
                    .filter_map(|r| r.cost.s)
                    .collect();
                let mean = stats::mean(&ok);
                // Strata: this arm's verified S per workload.
                let mut by_w: std::collections::BTreeMap<&str, Vec<f64>> =
                    std::collections::BTreeMap::new();
                for r in mine.iter().filter(|r| r.counts()) {
                    if let Some(v) = r.cost.s {
                        by_w.entry(r.workload.workload_id.as_str())
                            .or_default()
                            .push(v);
                    }
                }
                let strata: Vec<Vec<f64>> = by_w.into_values().collect();
                let ci = stats::cluster_bootstrap_ci(&strata, 4000, 20261002, |s| {
                    stats::mean(&s.concat())
                })
                .or_else(|| {
                    stats::bootstrap_ci(&strata, 4000, 20261002, |s| stats::mean(&s.concat()))
                });
                let any = mine.first();
                let floor = any.map(|r| r.boundaries.floor_s).unwrap_or(f64::NAN);
                let theory = any.and_then(|r| {
                    methods::expected_s(&arm.method.id, r.boundaries.automorphisms_available)
                });
                // S / reference over matched (workload, round) pairs, both
                // verified: a ratio of sums, as `compare` forms it.
                let ref_ratio = ref_arm.as_ref().and_then(|ra| {
                    let refs: std::collections::BTreeMap<(&str, u32), f64> = recs
                        .iter()
                        .filter(|r| r.counts() && &r.arm == ra && r.workload.curve.slug == *slug)
                        .filter_map(|r| {
                            Some(((r.workload.workload_id.as_str(), r.round), r.cost.s?))
                        })
                        .collect();
                    let (mut sa, mut sb) = (0.0, 0.0);
                    for r in mine.iter().filter(|r| r.counts()) {
                        if let (Some(v), Some(rv)) = (
                            r.cost.s,
                            refs.get(&(r.workload.workload_id.as_str(), r.round)),
                        ) {
                            sa += rv;
                            sb += v;
                        }
                    }
                    (sa > 0.0).then(|| sb / sa)
                });
                let mut levels: std::collections::BTreeMap<&str, u32> =
                    std::collections::BTreeMap::new();
                for r in &mine {
                    *levels.entry(r.isolation.level.name()).or_insert(0) += 1;
                }
                let o = |v: Option<f64>| v.map(|x| format!("{x:.3}")).unwrap_or_else(|| "–".into());
                println!(
                    "| {} | {} | {} | {} | {:.2} | {}/{} | {} | {} | {} | {} | {} | {} | {} | {} |",
                    s.session_id,
                    arm.name,
                    arm.method.id,
                    slug,
                    any.map(|r| (r.workload.curve.r as f64).log2())
                        .unwrap_or(f64::NAN),
                    ok.len(),
                    mine.len(),
                    o(mean),
                    ci.map(|(l, h)| format!("[{l:.3}, {h:.3}]"))
                        .unwrap_or_else(|| "–".into()),
                    o(theory),
                    o(match (mean, theory) {
                        (Some(m), Some(t)) => Some(m / t),
                        _ => None,
                    }),
                    o(mean.map(|m| m / floor)),
                    o(ref_ratio),
                    if mine.iter().any(|r| r.cost.lower_bound) {
                        "yes"
                    } else {
                        "no"
                    },
                    levels
                        .iter()
                        .map(|(k, v)| format!("{k}×{v}"))
                        .collect::<Vec<_>>()
                        .join(" "),
                );
            }
        }
    }
    Ok(())
}
