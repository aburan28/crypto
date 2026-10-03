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
//! ecbench db sql PATHS... | sqlite3 ecbench.db
//! ecbench isolab-job --spec S --binary B    an isolab.job/v1 for an independent runner
//! ```
//!
//! The standard is `docs/ecbench/README.md`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};

use crypto_lib::cryptanalysis::ecbench::{
    audit, compare, db, host, isolab, isolation, methods, record, runner, spec, stats, workload,
};

#[derive(Parser)]
#[command(
    name = "ecbench",
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
    /// Audit a session from its own files.
    Verify {
        #[arg(long)]
        dir: String,
        /// Re-execute this many measured runs and require identical counts.
        #[arg(long, default_value_t = 0)]
        replay: usize,
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

fn print_json(v: &impl serde::Serialize) -> Result<(), String> {
    println!(
        "{}",
        serde_json::to_string_pretty(v).map_err(|e| e.to_string())?
    );
    Ok(())
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
                    "commit      {} (dirty: {:?})",
                    c.build.git_commit.as_deref().unwrap_or("?"),
                    c.build.git_dirty
                );
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
                exe: std::env::current_exe().map_err(|e| e.to_string())?,
            };
            let s = runner::run_session(&text, &out, opts, |line| {
                if !quiet {
                    eprintln!("{line}");
                }
            })?;
            eprintln!(
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
        Cmd::Verify {
            dir,
            replay,
            out,
            exit_code,
        } => {
            let dir = PathBuf::from(isolab::expand_env(&dir)?);
            let rep = audit::audit(&dir, replay)?;
            let text = serde_json::to_string_pretty(&rep).map_err(|e| e.to_string())? + "\n";
            if let Some(o) = out {
                std::fs::write(&o, &text).map_err(|e| format!("{}: {e}", o.display()))?;
                eprintln!(
                    "receipt sha256 {}",
                    crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex(text.as_bytes())
                );
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
    }
    Ok(ExitCode::SUCCESS)
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
                let refm = ref_arm.as_ref().and_then(|ra| {
                    let v: Vec<f64> = recs
                        .iter()
                        .filter(|r| r.counts() && &r.arm == ra && r.workload.curve.slug == *slug)
                        .filter_map(|r| r.cost.s)
                        .collect();
                    stats::mean(&v)
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
                    o(match (mean, refm) {
                        (Some(m), Some(r)) if r > 0.0 => Some(m / r),
                        _ => None,
                    }),
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
