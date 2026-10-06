//! Native source-bound SAT target producer; a frozen controller must grant
//! the one-use claim and independently audit the retained execution.
#[path = "prepared_sat_target_worker/contract.rs"]
mod contract;
#[path = "prepared_target_worker/journal.rs"]
mod journal;
#[path = "prepared_sat_worker/native.rs"]
#[allow(dead_code)]
mod native;
#[path = "prepared_ordinary_worker/contract.rs"]
#[allow(dead_code)]
mod preparation;

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::{
    prepared_n17_target::{PreparedN17Target, SatTargetBackend, SatTargetPlan},
    prepared_sat_control::{canonical_sha, sha256, QueryOutput},
};
use serde_json::{json, Value};
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
    BuildIdentity,
    Run {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        card: PathBuf,
        #[arg(long)]
        registration_sha256: String,
    },
    /// Internal parent/PID handshake for an accepted native role only.
    Child {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        card: PathBuf,
        #[arg(long)]
        attempt: PathBuf,
        #[arg(long)]
        role: String,
    },
    /// Retain structure of the original attempt prefix; never admit runtime.
    Inspect {
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
        #[arg(long)]
        out: PathBuf,
    },
}
fn self_binding(capsule: &Path, record: &contract::Registration) -> Result<PathBuf, String> {
    contract::require_identity(&contract::identity(), &record.worker_build_identity)?;
    let own = std::env::current_exe().map_err(|e| e.to_string())?;
    native::require(
        own.canonicalize().map_err(|e| e.to_string())?
            == capsule
                .join("immutable/bin")
                .join(contract::WORKER)
                .canonicalize()
                .map_err(|e| e.to_string())?
            && sha256(&native::read(&own, 128 * 1024 * 1024)?) == record.worker_sha256,
        "running SAT target worker differs from frozen executable",
    )?;
    Ok(own)
}
fn recheck_preparation(
    record: &contract::Registration,
    cfg: &contract::Config,
    execution: &Path,
) -> Result<Value, String> {
    let (root, original, auditor) = contract::preparation_paths(&record.preparation)?;
    let out = execution.join("preparation-recheck.json");
    let args = vec![
        "ordinary-control-audit".into(),
        "--capsule".into(),
        root.to_string_lossy().into_owned(),
        "--execution".into(),
        original.to_string_lossy().into_owned(),
        "--registration-sha256".into(),
        record.preparation.registration_sha256.clone(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
    ];
    let environment = contract::environment().into_iter().collect::<Vec<_>>();
    let result = native::measured_child_request(native::ChildRequest {
        program: &auditor,
        args: &args,
        cwd: &root,
        stem: &execution.join("preparation-auditor"),
        deadline_ms: cfg.worker_timeout_ms,
        ledger: &execution.join("preparation-auditor-pids"),
        helper: false,
        input: None,
        environment: &environment,
    })?;
    native::require(
        result.exit_code == Some(0) && !result.timed_out,
        "SAT original preparation auditor failed; output retained",
    )?;
    let producer_bytes = native::read(&original.join("producer.json"), 16 * 1024 * 1024)?;
    let producer: Value = serde_json::from_slice(&producer_bytes).map_err(|e| e.to_string())?;
    let math = contract::preparation_receipt(
        &record.preparation,
        &native::load(&out)?,
        &producer,
        &producer_bytes,
    )?;
    contract::preparation_paths(&record.preparation)?;
    native::require(
        native::read(&original.join("producer.json"), 16 * 1024 * 1024)? == producer_bytes,
        "SAT preparation producer changed during recheck",
    )?;
    Ok(math)
}
struct Backend {
    capsule: PathBuf,
    execution: PathBuf,
    card: PathBuf,
    worker: PathBuf,
    cfg: contract::Config,
    journal: journal::Journal,
}
impl SatTargetBackend for Backend {
    fn started(&mut self, row: &Value) -> Result<(), String> {
        self.journal.start(row)
    }
    fn query(
        &mut self,
        trial: usize,
        point: [u64; 2],
        plan: &SatTargetPlan,
    ) -> Result<QueryOutput, String> {
        native::require(
            trial < self.cfg.plan.max_queries
                && point.iter().all(|&v| v < 1 << 17)
                && json!(plan) == json!(self.cfg.plan),
            "SAT target native query differs from frozen plan",
        )?;
        let dir = self.execution.join(format!("query-{trial:03}"));
        fs::create_dir(&dir).map_err(|e| e.to_string())?;
        native::save(
            &dir.join("query.json"),
            &json!(contract::Query { trial, point }),
        )?;
        contract::check_roles(&self.capsule)?;
        let args = |role: &str| {
            vec![
                "child".into(),
                "--capsule".into(),
                self.capsule.to_string_lossy().into_owned(),
                "--card".into(),
                self.card.to_string_lossy().into_owned(),
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
            plan.exporter_timeout_ms,
            &ledger,
            true,
        )?;
        native::require(
            exporter.exit_code == Some(0) && !exporter.timed_out,
            "SAT target exporter failed; original output retained",
        )?;
        let (manifest, anf, cnf, before) = native::validate_exports(&dir.join("instance"), point)?;
        let output = native::measured_child(
            &self.worker,
            &args("cms"),
            &dir,
            &dir.join("cms"),
            plan.solver_timeout_ms,
            &ledger,
            true,
        )?;
        let (manifest_after, anf_after, cnf_after, after) =
            native::validate_exports(&dir.join("instance"), point)?;
        native::require(
            manifest == manifest_after && anf == anf_after && cnf == cnf_after && before == after,
            "SAT target native source changed during solver run",
        )?;
        contract::check_roles(&self.capsule)?;
        Ok(QueryOutput {
            native: output,
            manifest,
            anf,
            cnf,
            source_receipt: json!({"files":before,"files_unchanged_before_after":true,
                "exporter":exporter.receipt}),
        })
    }
    fn completed(&mut self, row: &Value) -> Result<(), String> {
        self.journal.complete(row)
    }
}
fn run(capsule: &Path, execution: &Path, card_path: &Path, seal: &str) -> Result<(), String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let card_path = card_path.canonicalize().map_err(|e| e.to_string())?;
    let record = contract::check_capsule(&capsule, seal)?;
    let (card, card_sha256) = contract::claim(&capsule, &execution, &card_path, &record, seal)?;
    native::require(
        std::env::vars().collect::<std::collections::BTreeMap<_, _>>() == contract::environment(),
        "SAT target worker environment differs from frozen envelope",
    )?;
    let worker = self_binding(&capsule, &record)?;
    let cfg = contract::config(&capsule)?;
    native::save(
        &execution.join("worker-started.json"),
        &json!({"schema_version":1,"scope":contract::SCOPE,"registration_sha256":seal,
            "worker_sha256":record.worker_sha256,"worker_build_identity":contract::identity(),
            "config_sha256":record.config_sha256,"card_sha256":card_sha256,
            "target":card.target,"pid":std::process::id(),
            "source_bound_execution_admitted":false}),
    )?;
    let result = (|| -> Result<Value, String> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .map_err(|e| e.to_string())?;
        let math = recheck_preparation(&record, &cfg, &execution)?;
        let prepared = PreparedN17Target::from_ordinary_math(&math)?;
        let journal = journal::Journal::create(
            &execution,
            cfg.journal_binding(seal, &record.worker_sha256, card.target),
        )?;
        let plan = cfg.plan.clone();
        let mut backend = Backend {
            capsule: capsule.clone(),
            execution: execution.clone(),
            card: card_path.clone(),
            worker,
            cfg,
            journal,
        };
        prepared.solve_external_sat(card.target, &plan, &mut backend)
    })();
    let child_drain = native::drain_ledger(&execution.join("child-pids"));
    let prep_drain = native::drain_ledger(&execution.join("preparation-auditor-pids"));
    let unchanged = contract::check_capsule(&capsule, seal);
    let card_unchanged = contract::card(&card_path).map(|(_, hash)| hash == card_sha256);
    let recovered = result.as_ref().is_ok_and(|report| {
        report["verified_scalar"].is_u64()
            && report["failure"].is_null()
            && report["independent_recovery_certificate"]["inside_online_interval"] == true
    });
    match &result {
        Ok(value) => native::save(&execution.join("producer.json"), value)?,
        Err(reason) => native::save(
            &execution.join("producer-error.json"),
            &json!({"reason":reason,"source_bound_execution_admitted":false,
                "full_goal_complete":false}),
        )?,
    }
    native::save(
        &execution.join("worker-terminal.json"),
        &json!({"schema_version":1,"scope":contract::SCOPE,"registration_sha256":seal,
            "card_sha256":card_sha256,"child_drain_confirmed":child_drain.is_ok(),
            "child_drain_error":child_drain.as_ref().err(),
            "preparation_auditor_drain_confirmed":prep_drain.is_ok(),
            "preparation_auditor_drain_error":prep_drain.as_ref().err(),
            "immutable_binding_unchanged":unchanged.is_ok(),
            "immutable_binding_error":unchanged.as_ref().err(),
            "card_unchanged":card_unchanged.as_ref().is_ok_and(|v| *v),
            "card_error":card_unchanged.as_ref().err(),
            "producer_verified_recovery":recovered,
            "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
            "headline_eligible":false,"promotion_eligible":false,
            "full_goal_complete":false,"online_speedup":null}),
    )?;
    result?;
    child_drain?;
    prep_drain?;
    unchanged?;
    native::require(card_unchanged?, "SAT target card changed during execution")?;
    native::require(recovered, "SAT target did not recover a verified scalar")?;
    Ok(())
}
fn child(capsule: &Path, card: &Path, attempt: &Path, role: &str) -> Result<(), String> {
    native::enforce_hardware()?;
    let mut handshake = String::new();
    std::io::stdin()
        .take(4)
        .read_to_string(&mut handshake)
        .map_err(|e| e.to_string())?;
    native::require(
        handshake == "GO\n",
        "SAT target helper lacks launch handshake",
    )?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    #[cfg(unix)]
    let (cfg, query, claimed) =
        contract::helper_context(&capsule, card, attempt, unsafe { libc::getppid() } as u32)?;
    #[cfg(not(unix))]
    return Err("SAT target native helper requires Unix".into());
    #[cfg(unix)]
    {
        let record = contract::registration(&capsule, &claimed.registration_sha256)?;
        self_binding(&capsule, &record)?;
        native::require(
            unsafe { libc::setpgid(0, 0) } == 0,
            "SAT target helper process group isolation failed",
        )?;
        let (program, args, pin) = contract::role(&capsule, &cfg, &query, role)?;
        native::require(
            sha256(&native::read(&program, 4 * 1024 * 1024)?) == pin,
            "SAT target native role changed before exec",
        )?;
        native::save(
            &attempt.join(format!("{role}.launch.json")),
            &json!({"role":role,"program":program,"argv":args,
                "binary_sha256_before_exec":pin,"config_sha256":canonical_sha(&json!(cfg))?,
                "query_sha256":canonical_sha(&json!(query))?,
                "card_sha256":claimed.card_sha256,
                "environment":{"LC_ALL":"C"},"cwd":attempt,"pid":std::process::id()}),
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
            serde_json::to_string_pretty(&contract::identity()).map_err(|e| e.to_string())
        }
        Action::Run {
            capsule,
            execution,
            card,
            registration_sha256,
        } => run(&capsule, &execution, &card, &registration_sha256)
            .map(|_| "SAT target producer terminated; independent admission pending".into()),
        Action::Child {
            capsule,
            card,
            attempt,
            role,
        } => child(&capsule, &card, &attempt, &role).map(|_| String::new()),
        Action::Inspect {
            execution,
            registration_sha256,
            out,
        } => journal::inspect(&execution, &registration_sha256).and_then(|retained| {
            native::require(
                out.parent()
                    .ok_or("inspection output lacks parent")?
                    .canonicalize()
                    .map_err(|e| e.to_string())?
                    != execution
                        .join("attempts")
                        .canonicalize()
                        .map_err(|e| e.to_string())?,
                "SAT target inspection output would alter journal",
            )?;
            native::save(&out, &retained)?;
            serde_json::to_string_pretty(&retained).map_err(|e| e.to_string())
        }),
    };
    match result {
        Ok(output) => {
            if !output.is_empty() {
                println!("{output}");
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}
