//! Run, replay and optionally publish a coordination receipt for one bounded
//! P-256 isogeny-walk task.

use clap::{ArgAction, Parser, Subcommand};
use crypto_lib::cryptanalysis::p256_isogeny_task::{
    build_result, plan_taskq, verify_result, CairnWalkReceipt, P256IsogenyTask,
    P256IsogenyTaskResult, CERTIFICATE_FILE, RESULT_FILE, TASK_FILE,
};
use crypto_lib::cryptanalysis::p256_isogeny_walk::{
    from_json_lines, generate, to_json_lines, WalkCertificate, DEFAULT_STEPS,
};
use crypto_lib::cryptanalysis::pollard_collab::cairn::{
    wall_clock, CairnConfig, CairnTransport, Submitter,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use serde::Serialize;
use serde_json::json;
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(about = "Portable, fully replayed work unit for the bounded P-256 isogeny walk")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Write a validated immutable task manifest.
    Manifest {
        #[arg(long)]
        run_id: String,
        #[arg(long)]
        task_id: String,
        #[arg(long)]
        source_commit: String,
        #[arg(long, default_value_t = DEFAULT_STEPS)]
        steps: u64,
        #[arg(long)]
        output: PathBuf,
    },
    /// Generate, internally replay and write one task result.
    Run {
        #[arg(long, conflicts_with = "task_json", required_unless_present = "task_json")]
        task: Option<PathBuf>,
        #[arg(long, conflicts_with = "task", required_unless_present = "task")]
        task_json: Option<String>,
        /// Output directory; defaults to TASKQ_OUTPUT_DIR.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Independently replay all result files and compare every digest/count.
    Verify {
        #[arg(long)]
        task: PathBuf,
        #[arg(long)]
        dir: PathBuf,
    },
    /// Emit one taskq spec. This never submits it.
    PlanTaskq {
        #[arg(long)]
        task: PathBuf,
        #[arg(long, default_value = "crypto")]
        repo: String,
        #[arg(long, default_value = "cpu")]
        queue: String,
        #[arg(long, default_value_t = 1800)]
        timeout_seconds: u64,
        #[arg(long)]
        output: PathBuf,
    },
    /// Derive the deterministic, coordination-only Cairn receipt offline.
    CairnReceipt {
        #[arg(long)]
        task: PathBuf,
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        output: PathBuf,
    },
    /// Commit/reveal the coordination receipt. Run again after an epoch to reveal.
    PublishCairn {
        #[arg(long)]
        task: PathBuf,
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        node: String,
        #[arg(long)]
        objective: String,
        #[arg(long, conflicts_with = "submitter", required_unless_present = "submitter")]
        identity: Option<PathBuf>,
        #[arg(long, conflicts_with = "identity", required_unless_present = "identity")]
        submitter: Option<String>,
        #[arg(long)]
        state: PathBuf,
        #[arg(long, default_value_t = 600)]
        epoch_seconds: u64,
        /// Required acknowledgement: the receipt is not a payable-work proof.
        #[arg(long, action = ArgAction::SetTrue)]
        ack_coordination_only: bool,
    },
}

struct LoadedResult {
    result: P256IsogenyTaskResult,
    certificate: WalkCertificate,
    artifact_bytes: Vec<u8>,
}

fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|error| format!("parse {}: {error}", path.display()))
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if let Some(parent) = path.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| format!("create-new {}: {error}", path.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write {}: {error}", path.display()))
}

fn write_json_new(path: &Path, value: &impl Serialize) -> Result<(), String> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    bytes.push(b'\n');
    write_new(path, &bytes)
}

fn compress(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .write(Vec::new(), Compression::best());
    encoder
        .write_all(bytes)
        .map_err(|error| format!("compress certificate: {error}"))?;
    encoder
        .finish()
        .map_err(|error| format!("finish certificate: {error}"))
}

fn decompress(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut output = Vec::new();
    GzDecoder::new(bytes)
        .read_to_end(&mut output)
        .map_err(|error| format!("decompress certificate: {error}"))?;
    Ok(output)
}

fn task_from_inputs(
    task: Option<PathBuf>,
    task_json: Option<String>,
) -> Result<P256IsogenyTask, String> {
    let task = match (task, task_json) {
        (Some(path), None) => read_json(&path)?,
        (None, Some(text)) => serde_json::from_str(&text)
            .map_err(|error| format!("parse --task-json: {error}"))?,
        _ => return Err("pass exactly one of --task or --task-json".into()),
    };
    task.validate()?;
    Ok(task)
}

fn result_dir(out: Option<PathBuf>) -> Result<PathBuf, String> {
    out.or_else(|| std::env::var_os("TASKQ_OUTPUT_DIR").map(PathBuf::from))
        .ok_or_else(|| "--out or TASKQ_OUTPUT_DIR is required".into())
}

fn preflight_fresh(dir: &Path) -> Result<(), String> {
    fs::create_dir_all(dir).map_err(|error| format!("create {}: {error}", dir.display()))?;
    for name in [TASK_FILE, CERTIFICATE_FILE, RESULT_FILE] {
        let path = dir.join(name);
        if path.exists() {
            return Err(format!("refusing to overwrite {}", path.display()));
        }
    }
    Ok(())
}

fn load_verified(task: &P256IsogenyTask, dir: &Path) -> Result<LoadedResult, String> {
    task.validate()?;
    let stored_task: P256IsogenyTask = read_json(&dir.join(TASK_FILE))?;
    if &stored_task != task {
        return Err("result task.json does not match the requested task".into());
    }
    let result: P256IsogenyTaskResult = read_json(&dir.join(RESULT_FILE))?;
    let artifact_bytes = fs::read(dir.join(CERTIFICATE_FILE))
        .map_err(|error| format!("read {}/{}: {error}", dir.display(), CERTIFICATE_FILE))?;
    let certificate = from_json_lines(&decompress(&artifact_bytes)?)?;
    verify_result(task, &result, &artifact_bytes, &certificate)?;
    Ok(LoadedResult {
        result,
        certificate,
        artifact_bytes,
    })
}

fn receipt_for(task: &P256IsogenyTask, dir: &Path) -> Result<CairnWalkReceipt, String> {
    let loaded = load_verified(task, dir)?;
    CairnWalkReceipt::from_verified(
        task,
        &loaded.result,
        &loaded.artifact_bytes,
        &loaded.certificate,
    )
}

fn run() -> Result<(), String> {
    match Cli::parse().command {
        Command::Manifest {
            run_id,
            task_id,
            source_commit,
            steps,
            output,
        } => {
            let task = P256IsogenyTask::new(run_id, task_id, source_commit, steps);
            task.validate()?;
            write_json_new(&output, &task)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&task).map_err(|error| error.to_string())?
            );
        }
        Command::Run {
            task,
            task_json,
            out,
        } => {
            let task = task_from_inputs(task, task_json)?;
            let out = result_dir(out)?;
            preflight_fresh(&out)?;

            let certificate = generate(task.steps)?;
            let artifact_bytes = compress(&to_json_lines(&certificate)?)?;
            let result = build_result(&task, &artifact_bytes, &certificate)?;

            write_json_new(&out.join(TASK_FILE), &task)?;
            write_new(&out.join(CERTIFICATE_FILE), &artifact_bytes)?;
            write_json_new(&out.join(RESULT_FILE), &result)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&result).map_err(|error| error.to_string())?
            );
        }
        Command::Verify { task, dir } => {
            let task: P256IsogenyTask = read_json(&task)?;
            let loaded = load_verified(&task, &dir)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&loaded.result)
                    .map_err(|error| error.to_string())?
            );
        }
        Command::PlanTaskq {
            task,
            repo,
            queue,
            timeout_seconds,
            output,
        } => {
            let task: P256IsogenyTask = read_json(&task)?;
            let spec = plan_taskq(&task, &repo, &queue, timeout_seconds)?;
            write_json_new(&output, &spec)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&spec).map_err(|error| error.to_string())?
            );
        }
        Command::CairnReceipt { task, dir, output } => {
            let task: P256IsogenyTask = read_json(&task)?;
            let receipt = receipt_for(&task, &dir)?;
            write_json_new(&output, &receipt)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
            );
        }
        Command::PublishCairn {
            task,
            dir,
            node,
            objective,
            identity,
            submitter,
            state,
            epoch_seconds,
            ack_coordination_only,
        } => {
            if !ack_coordination_only {
                return Err(
                    "--ack-coordination-only is required; this receipt is not payable-work proof"
                        .into(),
                );
            }
            let task: P256IsogenyTask = read_json(&task)?;
            let receipt = receipt_for(&task, &dir)?;
            let who = match (identity, submitter) {
                (Some(path), None) => Submitter::from_identity_file(&path)?,
                (None, Some(name)) if !name.trim().is_empty() => Submitter::Nickname(name),
                _ => return Err("pass exactly one nonempty --identity or --submitter".into()),
            };
            if let Some(parent) = state.parent().filter(|parent| !parent.as_os_str().is_empty()) {
                fs::create_dir_all(parent)
                    .map_err(|error| format!("create {}: {error}", parent.display()))?;
            }
            let mut transport = CairnTransport::open(CairnConfig {
                url: node.trim_end_matches('/').to_string(),
                objective_id: objective,
                submitter: who,
                answer_objective: None,
                epoch_secs: epoch_seconds,
                state_path: Some(state),
                clock: wall_clock(),
            })?;
            let (revealed, refused) = transport.reveal_pending()?;
            let published = transport.publish_artifact(receipt.artifact()?, receipt.claim_key())?;
            println!(
                "{}",
                serde_json::to_string_pretty(&json!({
                    "schema": "p256.isogeny.cairn.publish-report/v1",
                    "committed": published.committed,
                    "skipped": published.skipped,
                    "revealed": revealed,
                    "refused": refused,
                    "pending": transport.pending(),
                    "submitter": transport.submitter(),
                    "work_sha256": receipt.work_sha256,
                    "claim_class": receipt.claim_class,
                }))
                .map_err(|error| error.to_string())?
            );
        }
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p256-isogeny-task: {error}");
        std::process::exit(1);
    }
}
