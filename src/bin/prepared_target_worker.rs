//! Frozen single public-point F5 worker, bounded to the disclosed n17 fixture.
#[path = "prepared_target_worker/contract.rs"]
mod contract;
#[path = "prepared_target_worker/journal.rs"]
mod journal;
#[path = "prepared_sat_worker/native.rs"]
#[allow(dead_code)] // Reuse transport/files only; never the old capsule or helper gate.
mod native;
#[path = "prepared_ordinary_worker/contract.rs"]
#[allow(dead_code)] // Verify the original preparation; never launch its worker.
mod preparation;
use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::{
    prepared_n17_target::PreparedN17Target, prepared_sat_control::sha256,
};
use serde_json::{json, Value};
use std::{
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
    /// Report compiled source/build identity; no preparation or solver.
    BuildIdentity,
    /// Requires a separate sealed capsule and already-consumed one-use claim.
    Run {
        #[arg(long)]
        capsule: PathBuf,
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
    },
    /// Read bounded attempt files as data; never runtime or mathematical admission.
    Inspect {
        #[arg(long)]
        execution: PathBuf,
        #[arg(long)]
        registration_sha256: String,
        #[arg(long)]
        out: PathBuf,
    },
}
fn self_binding(capsule: &Path, record: &contract::Registration) -> Result<(), String> {
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
        "target running executable differs from frozen worker",
    )
}
fn recheck_preparation(
    record: &contract::Registration,
    cfg: &contract::Config,
    execution: &Path,
) -> Result<Value, String> {
    let (capsule, original, checker) = contract::preparation_paths(&record.preparation)?;
    let out = execution.join("preparation-recheck.json");
    let args = vec![
        "ordinary-control-audit".into(),
        "--capsule".into(),
        capsule.to_string_lossy().into_owned(),
        "--execution".into(),
        original.to_string_lossy().into_owned(),
        "--registration-sha256".into(),
        record.preparation.registration_sha256.clone(),
        "--out".into(),
        out.to_string_lossy().into_owned(),
    ];
    let environment = contract::environment().into_iter().collect::<Vec<_>>();
    let result = native::measured_child_request(native::ChildRequest {
        program: &checker,
        args: &args,
        cwd: &capsule,
        stem: &execution.join("preparation-auditor"),
        deadline_ms: cfg.worker_timeout_ms,
        ledger: &execution.join("preparation-auditor-pids"),
        helper: false,
        input: None,
        environment: &environment,
    })?;
    native::require(
        result.exit_code == Some(0) && !result.timed_out,
        "original preparation audit failed; original output retained",
    )?;
    let producer_bytes = native::read(&original.join("producer.json"), 16 * 1024 * 1024)?;
    let producer: Value = serde_json::from_slice(&producer_bytes).map_err(|e| e.to_string())?;
    let math = contract::preparation_receipt(
        &record.preparation,
        &native::load(&out)?,
        &producer,
        &producer_bytes,
    )?;
    // Recheck after the auditor and original data read; no copied archive is run.
    contract::preparation_paths(&record.preparation)?;
    native::require(
        native::read(&original.join("producer.json"), 16 * 1024 * 1024)? == producer_bytes,
        "original preparation producer changed during recheck",
    )?;
    Ok(math)
}
fn run(capsule: &Path, execution: &Path, expected: &str) -> Result<(), String> {
    native::enforce_hardware()?;
    let capsule = capsule.canonicalize().map_err(|e| e.to_string())?;
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let record = contract::check_capsule(&capsule, expected)?;
    contract::claim(&capsule, &execution, &record, expected)?;
    native::require(
        std::env::vars().collect::<std::collections::BTreeMap<_, _>>() == contract::environment(),
        "target worker environment differs from frozen envelope",
    )?;
    self_binding(&capsule, &record)?;
    let cfg = contract::config(&capsule)?;
    native::save(
        &execution.join("worker-started.json"),
        &json!({"schema_version":1,"scope":contract::SCOPE,
        "registration_sha256":expected,"worker_sha256":record.worker_sha256,"worker_build_identity":contract::identity(),
        "config_sha256":record.config_sha256,"pid":std::process::id(),"source_bound_execution_admitted":false}),
    )?;
    let result = (|| -> Result<Value, String> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()
            .map_err(|e| e.to_string())?;
        let math = recheck_preparation(&record, &cfg, &execution)?;
        let prepared = PreparedN17Target::from_ordinary_math(&math)?;
        let mut journal =
            journal::Journal::create(&execution, cfg.binding(expected, &record.worker_sha256))?;
        prepared.solve_f5_recorded(cfg.target, &cfg.plan, &mut journal)
    })();
    let drain = native::drain_ledger(&execution.join("preparation-auditor-pids"));
    let unchanged = contract::check_capsule(&capsule, expected);
    let producer_verified_recovery = result.as_ref().is_ok_and(|r| {
        r["verified_scalar"].is_u64()
            && r["failure"].is_null()
            && r["independent_recovery_certificate"]["inside_online_interval"] == true
    });
    // A producer answer is provisional until the separately frozen target audit.
    match &result {
        Ok(value) => native::save(&execution.join("producer.json"), value)?,
        Err(reason) => native::save(
            &execution.join("producer-error.json"),
            &json!({"reason":reason,
            "source_bound_execution_admitted":false,"full_goal_complete":false}),
        )?,
    }
    native::save(
        &execution.join("worker-terminal.json"),
        &json!({"schema_version":1,"scope":contract::SCOPE,
        "registration_sha256":expected,"preparation_auditor_drain_confirmed":drain.is_ok(),
        "preparation_auditor_drain_error":drain.as_ref().err(),"immutable_binding_unchanged":unchanged.is_ok(),
        "immutable_binding_error":unchanged.as_ref().err(),"source_bound_execution_admitted":false,
        "producer_verified_recovery":producer_verified_recovery,
        "fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_speedup":null}),
    )?;
    result?;
    drain?;
    unchanged?;
    native::require(
        producer_verified_recovery,
        "target producer did not recover a verified scalar; original report retained",
    )?;
    Ok(())
}
fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Action::BuildIdentity => {
            serde_json::to_string_pretty(&contract::identity()).map_err(|e| e.to_string())
        }
        Action::Run {
            capsule,
            execution,
            registration_sha256,
        } => run(&capsule, &execution, &registration_sha256)
            .map(|_| "target producer terminated; independent admission pending".into()),
        Action::Inspect {
            execution,
            registration_sha256,
            out,
        } => journal::inspect(&execution, &registration_sha256).and_then(|v| {
            native::require(
                out.parent()
                    .ok_or("inspection output lacks parent")?
                    .canonicalize()
                    .map_err(|e| e.to_string())?
                    != execution
                        .join("attempts")
                        .canonicalize()
                        .map_err(|e| e.to_string())?,
                "inspection output would alter the attempt journal",
            )?;
            native::save(&out, &v)?;
            serde_json::to_string_pretty(&v).map_err(|e| e.to_string())
        }),
    };
    match result {
        Ok(output) => {
            println!("{output}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}
