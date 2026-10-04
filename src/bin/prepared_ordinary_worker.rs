//! Frozen target-free preparation worker for the disclosed educational n17 curve.
//! No registration/controller is implied by the API or its unit controls.
#[path = "prepared_ordinary_worker/contract.rs"]
mod contract;
#[path = "prepared_sat_worker/native.rs"]
#[allow(dead_code)] // Reuse only the bounded transport, not its old capsule/helper gate.
mod native;
use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::{
    koblitz_index_calculus::probe_scalar,
    prepared_ordinary::{self as producer, Family, SatBackend},
    prepared_sat_control::{canonical_sha, sha256, QueryOutput},
};
use serde_json::json;
#[cfg(unix)]
use std::os::unix::process::CommandExt;
use std::{
    fs,
    io::Read,
    path::{Path, PathBuf},
    process::{Command, ExitCode, Stdio},
};

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Action,
}
#[derive(Subcommand)]
enum Action {
    /// Return compiled source/build identity; no preparation or solver.
    BuildIdentity,
    /// Requires a new matching frozen capsule and durable one-use claim.
    Run {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
    },
    /// Internal parent/PID handshake for the accepted native roles only.
    Child {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        attempt: PathBuf,
        #[arg(long)]
        role: String,
    },
}
struct Backend {
    capsule: PathBuf,
    execution: PathBuf,
    worker: PathBuf,
    cfg: contract::Config,
}
impl SatBackend for Backend {
    fn query(&mut self, trial: usize, point: [u64; 2]) -> Result<QueryOutput, String> {
        native::require(
            trial < self.cfg.plan.planned_queries,
            "ordinary backend exceeded frozen panel",
        )?;
        let dir = self.execution.join(format!("query-{trial:03}"));
        fs::create_dir(&dir).map_err(|e| e.to_string())?;
        native::save(
            &dir.join("query.json"),
            &json!(contract::Query {
                trial,
                scalar: probe_scalar(self.cfg.plan.algorithm_seed, trial as u64, 65587),
                point,
            }),
        )?;
        contract::check_roles(&self.capsule)?;
        let args = |role: &str| {
            vec![
                "child".into(),
                "--capsule".into(),
                self.capsule.to_string_lossy().into_owned(),
                "--attempt".into(),
                dir.to_string_lossy().into_owned(),
                "--role".into(),
                role.into(),
            ]
        };
        let ledger = self.execution.join("child-pids");
        let exporter = native::measured_child(
            &self.worker,
            &args("exporter"),
            &dir,
            &dir.join("exporter"),
            self.cfg.exporter_timeout_ms,
            &ledger,
            true,
        )?;
        native::require(
            exporter.exit_code == Some(0) && !exporter.timed_out,
            "ordinary exporter failed; original output retained",
        )?;
        let (manifest, anf, cnf, before) = native::validate_exports(&dir.join("instance"), point)?;
        let output = native::measured_child(
            &self.worker,
            &args("cms"),
            &dir,
            &dir.join("cms"),
            self.cfg.solver_timeout_ms,
            &ledger,
            true,
        )?;
        let (manifest_after, anf_after, cnf_after, after) =
            native::validate_exports(&dir.join("instance"), point)?;
        native::require(
            manifest == manifest_after && anf == anf_after && cnf == cnf_after && before == after,
            "ordinary native source changed during solve; original files retained",
        )?;
        contract::check_roles(&self.capsule)?;
        Ok(QueryOutput {
            native: output,
            manifest,
            anf,
            cnf,
            source_receipt: json!({"files":before,"files_unchanged_before_after":true,"exporter":exporter.receipt}),
        })
    }
}
fn self_binding(capsule: &Path, record: &contract::Registration) -> Result<PathBuf, String> {
    contract::require_identity(&contract::identity(), &record.worker_build_identity)?;
    let worker = std::env::current_exe().map_err(|e| e.to_string())?;
    native::require(
        worker.canonicalize().map_err(|e| e.to_string())?
            == capsule
                .join("immutable/bin")
                .join(contract::WORKER)
                .canonicalize()
                .map_err(|e| e.to_string())?
            && sha256(&native::read(&worker, 128 * 1024 * 1024)?) == record.worker_sha256,
        "ordinary running worker differs from frozen executable",
    )?;
    Ok(worker)
}
fn run(capsule: &Path, execution: &Path) -> Result<(), String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let (record, hash) = contract::check_capsule(&capsule)?;
    contract::claim(&capsule, &execution, &record, &hash)?;
    native::require(
        std::env::vars().collect::<std::collections::BTreeMap<_, _>>() == contract::environment(),
        "ordinary worker environment differs from frozen envelope",
    )?;
    let worker = self_binding(&capsule, &record)?;
    native::save(
        &execution.join("worker-started.json"),
        &json!({"registration_sha256":hash,
        "worker_sha256":record.worker_sha256,"worker_build_identity":contract::identity(),"pid":std::process::id()}),
    )?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .map_err(|e| e.to_string())?;
    let cfg = contract::config(&capsule)?;
    let mut observer = contract::DurableObserver::new(&execution);
    let report = match cfg.plan.family {
        Family::MatrixF5 => producer::prepare_f5(&cfg.plan, &mut observer),
        Family::Cryptominisat => {
            let mut backend = Backend {
                capsule: capsule.clone(),
                execution: execution.clone(),
                worker,
                cfg: cfg.clone(),
            };
            producer::prepare_sat(&cfg.plan, &mut backend, &mut observer)
        }
    };
    // Drain and immutable check still run after a producer error; neither
    // failure can turn a partial prefix into a scientific success.
    let drain = native::drain_ledger(&execution.join("child-pids"));
    let unchanged = contract::check_capsule(&capsule);
    match &report {
        Ok(report) => {
            native::save(
                &execution.join("mathematical-input.json"),
                &report["mathematical_input"],
            )?;
            native::save(&execution.join("producer.json"), report)?;
        }
        Err(reason) => native::save(
            &execution.join("producer-error.json"),
            &json!({"reason":reason,
            "source_bound_execution_admitted":false,"full_goal_complete":false}),
        )?,
    }
    native::save(
        &execution.join("worker-terminal.json"),
        &json!({
        "child_drain_confirmed":drain.is_ok(),"child_drain_error":drain.as_ref().err(),
        "immutable_binding_unchanged":unchanged.is_ok(),"immutable_binding_error":unchanged.as_ref().err(),
        "source_bound_execution_admitted":false,"full_goal_complete":false}),
    )?;
    report?;
    drain?;
    unchanged?;
    Ok(())
}
fn helper(capsule: &Path, attempt: &Path, role: &str) -> Result<(), String> {
    native::enforce_hardware()?;
    let mut handshake = String::new();
    std::io::stdin()
        .take(4)
        .read_to_string(&mut handshake)
        .map_err(|e| e.to_string())?;
    native::require(
        handshake == "GO\n",
        "ordinary helper lacks durable launch handshake",
    )?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    #[cfg(unix)]
    let (cfg, q) = contract::helper_context(&capsule, attempt, unsafe { libc::getppid() } as u32)?;
    #[cfg(not(unix))]
    return Err("ordinary native helper requires Unix".into());
    #[cfg(unix)]
    {
        let (record, _) = contract::registration(&capsule)?;
        self_binding(&capsule, &record)?;
        native::require(
            unsafe { libc::setpgid(0, 0) } == 0,
            "ordinary helper group isolation failed",
        )?;
        let (program, args, pin) = contract::role(&capsule, &cfg, &q, role)?;
        native::require(
            sha256(&native::read(&program, 4 * 1024 * 1024)?) == pin,
            "ordinary native role changed before exec",
        )?;
        native::save(
            &attempt.join(format!("{role}.launch.json")),
            &json!({"role":role,"program":program,"argv":args,
            "binary_sha256_before_exec":pin,"config_sha256":canonical_sha(&json!(cfg))?,
            "query_sha256":canonical_sha(&json!(q))?,"environment":{"LC_ALL":"C"},"cwd":attempt,"pid":std::process::id()}),
        )?;
        let mut command = Command::new(program);
        command
            .args(args)
            .current_dir(attempt)
            .env_clear()
            .env("LC_ALL", "C")
            .stdin(Stdio::null());
        Err(command.exec().to_string())
    }
}
fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Action::BuildIdentity => {
            println!("{}", contract::identity());
            Ok(())
        }
        Action::Run { capsule, execution } => run(&capsule, &execution),
        Action::Child {
            capsule,
            attempt,
            role,
        } => helper(&capsule, &attempt, &role),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
