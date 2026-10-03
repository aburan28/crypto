//! Registered native producer for one disclosed synthetic n17 control.
#[path = "prepared_sat_worker/native.rs"]
#[allow(dead_code)] // The same transport also serves registration and audit in icprog.
mod native;
use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::prepared_sat_control::{
    self as producer, ControlConfig, OnlineClock, QueryBackend, QueryOutput,
};
use serde_json::{json, Value};
use std::{
    fs,
    path::{Path, PathBuf},
    process::ExitCode,
};

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Action,
}
#[derive(Subcommand)]
enum Action {
    /// A one-use execution claim must already exist. No ambient executable choices.
    Run {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
    },
    /// Internal handshake helper. Accepts only the two pinned native roles.
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
}
impl QueryBackend for Backend {
    fn query(
        &mut self,
        trial: usize,
        point: [u64; 2],
        cfg: &ControlConfig,
        _clock: &mut OnlineClock,
    ) -> Result<QueryOutput, String> {
        let dir = self.execution.join(format!("query-{trial:02}"));
        fs::create_dir(&dir).map_err(|e| e.to_string())?;
        native::save(
            &dir.join("query.json"),
            &json!({"trial":trial,"point":point}),
        )?;
        let ledger = self.execution.join("child-pids");
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
        let exporter = native::measured_child(
            &self.worker,
            &args("exporter"),
            &dir,
            &dir.join("exporter"),
            cfg.exporter_timeout_ms,
            &ledger,
            true,
        )?;
        native::require(
            exporter.exit_code == Some(0) && !exporter.timed_out,
            "exporter failed; original output retained",
        )?;
        let (manifest, anf, cnf, source_receipt) =
            native::validate_exports(&dir.join("instance"), point)?;
        let native = native::measured_child(
            &self.worker,
            &args("cms"),
            &dir,
            &dir.join("cms"),
            cfg.solver_timeout_ms,
            &ledger,
            true,
        )?;
        let (_, anf_after, cnf_after, after) =
            native::validate_exports(&dir.join("instance"), point)?;
        native::require(
            after == source_receipt && anf == anf_after && cnf == cnf_after,
            "solver input changed during execution",
        )?;
        // Helpers exec the accepted role. Check role binaries separately from the helper.
        native::require(
            producer::sha256(&native::read(
                &self.capsule.join("immutable/assets/bin/cms"),
                4 * 1024 * 1024,
            )?) == native::CMS_SHA
                && producer::sha256(&native::read(
                    &self.capsule.join("immutable/assets/bin/exporter"),
                    2 * 1024 * 1024,
                )?) == native::EXPORTER_SHA,
            "native role binary changed during query",
        )?;
        Ok(QueryOutput {
            native,
            manifest,
            anf,
            cnf,
            source_receipt: json!({"files":source_receipt,"files_unchanged_before_after":true,"exporter":exporter.receipt}),
        })
    }
    fn progress(&mut self, attempt: &Value) -> Result<(), String> {
        let trial = attempt["trial"].as_u64().ok_or("missing progress trial")?;
        native::save(
            &self.execution.join(format!("progress-{trial:02}.json")),
            attempt,
        )
    }
}
fn run(capsule: &Path, execution: &Path) -> Result<(), String> {
    native::enforce_hardware()?;
    let registration = native::check_capsule(capsule)?;
    let claim = native::load(&capsule.join("consumed.json"))?;
    native::require(
        claim["registration_sha256"] == producer::canonical_sha(&registration)?
            && claim["execution"] == json!(execution),
        "one-use claim differs",
    )?;
    native::save(
        &execution.join("worker-started.json"),
        &json!({"registration_sha256":producer::canonical_sha(&registration)?,"pid":std::process::id()}),
    )?;
    let worker = std::env::current_exe().map_err(|e| e.to_string())?;
    native::require(
        producer::sha256(&native::read(&worker, 128 * 1024 * 1024)?)
            == registration["worker_sha256"]
                .as_str()
                .ok_or("worker pin missing")?,
        "running producer differs from capsule",
    )?;
    let prep = native::load(&capsule.join("preparation.json"))?;
    let cfg = native::config(capsule)?;
    let mut backend = Backend {
        capsule: capsule.into(),
        execution: execution.into(),
        worker,
    };
    let report = producer::solve(&prep, &cfg, &mut backend)?;
    native::drain_ledger(&execution.join("child-pids"))?;
    native::check_capsule(capsule)?;
    native::save(&execution.join("producer.json"), &report)?;
    Ok(())
}
fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Action::Run { capsule, execution } => run(&capsule, &execution),
        Action::Child {
            capsule,
            attempt,
            role,
        } => native::helper(&capsule, &attempt, &role),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
