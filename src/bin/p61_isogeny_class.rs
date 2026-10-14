//! Generate and independently replay bounded-memory p61 class segments.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, ErrorKind, Write};
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::million_store::{digest_file, FileDigest};
use crypto_lib::cryptanalysis::isogeny_walk::p61_class::{
    generate_segment_jsonl, plan_record, verify_segment_jsonl, Checkpoint, SegmentConfig,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use serde::{de::DeserializeOwned, Serialize};
use serde_json::json;

#[derive(Debug, Serialize)]
struct DurableWrite {
    digest: FileDigest,
    sync_status: String,
}

#[derive(Debug, Parser)]
#[command(about = "Generate and replay checkpointed p61 complete-class segments")]
struct Cli {
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Print the frozen class-group lattice and completion gates.
    Plan,
    /// Generate one immutable segment and its terminal checkpoint.
    GenerateSegment {
        #[arg(long)]
        x_start: u64,
        #[arg(long)]
        length: u64,
        #[arg(long, default_value_t = 256)]
        batch_steps: usize,
        #[arg(long, default_value_t = 1)]
        audit_points: usize,
        #[arg(long, default_value_t = 0)]
        audit_seed_x: u64,
        #[arg(long)]
        source_commit: String,
        #[arg(long)]
        class_number_receipt_sha256: String,
        #[arg(long)]
        previous_checkpoint: Option<PathBuf>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        checkpoint_output: PathBuf,
        #[arg(long)]
        receipt_output: PathBuf,
    },
    /// Replay every record and kernel in one segment.
    VerifySegment {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        checkpoint: PathBuf,
        #[arg(long)]
        previous_checkpoint: Option<PathBuf>,
        #[arg(long, default_value_t = 2)]
        audit_points: usize,
        #[arg(long, default_value_t = 7)]
        audit_seed_x: u64,
        #[arg(long)]
        receipt_output: PathBuf,
    },
}

fn partial_path(path: &Path) -> PathBuf {
    let mut name: OsString = path.as_os_str().to_owned();
    name.push(".partial");
    PathBuf::from(name)
}

fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|error| format!("parse {}: {error}", path.display()))
}

fn sync_file(file: &File, path: &Path) -> Result<String, String> {
    match file.sync_all() {
        Ok(()) => Ok("sync_all_completed".into()),
        Err(error)
            if error.kind() == ErrorKind::Unsupported || error.raw_os_error() == Some(45) =>
        {
            Ok("sync_all_unsupported_by_filesystem".into())
        }
        Err(error) => Err(format!("sync {}: {error}", path.display())),
    }
}

fn write_json_new<T: Serialize>(path: &Path, value: &T) -> Result<DurableWrite, String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    if let Some(parent) = path.parent().filter(|value| !value.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let partial = partial_path(path);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, value).map_err(|error| error.to_string())?;
    writer.write_all(b"\n").map_err(|error| error.to_string())?;
    writer.flush().map_err(|error| error.to_string())?;
    let sync_status = sync_file(writer.get_ref(), &partial)?;
    fs::rename(&partial, path).map_err(|error| {
        format!(
            "commit {} -> {}: {error}",
            partial.display(),
            path.display()
        )
    })?;
    Ok(DurableWrite {
        digest: digest_file(path)?,
        sync_status,
    })
}

fn checkpoint_with_digest(path: &Path) -> Result<(Checkpoint, FileDigest), String> {
    let checkpoint: Checkpoint = read_json(path)?;
    let digest = digest_file(path)?;
    Ok((checkpoint, digest))
}

#[allow(clippy::too_many_arguments)]
fn generate(
    x_start: u64,
    length: u64,
    batch_steps: usize,
    audit_points: usize,
    audit_seed_x: u64,
    source_commit: String,
    class_number_receipt_sha256: String,
    previous_checkpoint_path: Option<&Path>,
    output: &Path,
    checkpoint_output: &Path,
    receipt_output: &Path,
) -> Result<(), String> {
    if output.exists() || checkpoint_output.exists() || receipt_output.exists() {
        return Err("refusing to overwrite a segment, checkpoint, or receipt".into());
    }
    let previous = previous_checkpoint_path
        .map(checkpoint_with_digest)
        .transpose()?;
    let previous_digest = previous.as_ref().map(|(_, digest)| digest.sha256.clone());
    let config = SegmentConfig {
        x_start,
        length,
        batch_steps,
        audit_points,
        audit_seed_x,
        source_commit,
        class_number_receipt_sha256,
        previous_checkpoint_sha256: previous_digest,
    };
    config.validate()?;
    if let Some(parent) = output
        .parent()
        .filter(|value| !value.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let partial = partial_path(output);
    let file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let buffered = BufWriter::new(file);
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .operating_system(255)
        .write(buffered, Compression::new(6));
    let result = generate_segment_jsonl(
        &mut encoder,
        &config,
        previous
            .as_ref()
            .map(|(checkpoint, digest)| (checkpoint, digest.sha256.as_str())),
    );
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", partial.display()))?;
    buffered.flush().map_err(|error| error.to_string())?;
    let segment_sync_status = sync_file(buffered.get_ref(), &partial)?;
    let mut receipt = match result {
        Ok(value) => value,
        Err(error) => {
            return Err(format!(
                "{error}; preserved incomplete segment at {}",
                partial.display()
            ))
        }
    };
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    let segment_digest = digest_file(output)?;
    receipt.checkpoint.segment_sha256 = segment_digest.sha256.clone();
    let checkpoint_write = write_json_new(checkpoint_output, &receipt.checkpoint)?;
    let durable_receipt = json!({
        "status": "pass",
        "segment": output,
        "segment_digest": segment_digest,
        "checkpoint": checkpoint_output,
        "checkpoint_digest": checkpoint_write.digest,
        "durability": {
            "segment_sync_status": segment_sync_status,
            "checkpoint_sync_status": checkpoint_write.sync_status,
        },
        "receipt": receipt,
    });
    let receipt_write = write_json_new(receipt_output, &durable_receipt)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "receipt_output": receipt_output,
            "receipt_digest": receipt_write.digest,
            "receipt_sync_status": receipt_write.sync_status,
            "result": durable_receipt,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn verify(
    input: &Path,
    checkpoint_path: &Path,
    previous_checkpoint_path: Option<&Path>,
    audit_points: usize,
    audit_seed_x: u64,
    receipt_output: &Path,
) -> Result<(), String> {
    if receipt_output.exists() {
        return Err(format!(
            "refusing to overwrite {}",
            receipt_output.display()
        ));
    }
    let (checkpoint, checkpoint_digest) = checkpoint_with_digest(checkpoint_path)?;
    let segment_digest = digest_file(input)?;
    if checkpoint.segment_sha256 != segment_digest.sha256 {
        return Err("checkpoint does not bind the supplied segment".into());
    }
    let previous = previous_checkpoint_path
        .map(checkpoint_with_digest)
        .transpose()?;
    let file = File::open(input).map_err(|error| format!("open {}: {error}", input.display()))?;
    let reader = BufReader::new(GzDecoder::new(file));
    let receipt = verify_segment_jsonl(
        reader,
        &checkpoint,
        previous
            .as_ref()
            .map(|(prior, digest)| (prior, digest.sha256.as_str())),
        audit_points,
        audit_seed_x,
    )?;
    let durable_receipt = json!({
        "status": "pass",
        "segment": input,
        "segment_digest": segment_digest,
        "checkpoint": checkpoint_path,
        "checkpoint_digest": checkpoint_digest,
        "receipt": receipt,
    });
    let receipt_write = write_json_new(receipt_output, &durable_receipt)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "receipt_output": receipt_output,
            "receipt_digest": receipt_write.digest,
            "receipt_sync_status": receipt_write.sync_status,
            "result": durable_receipt,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn run() -> Result<(), String> {
    let cli = Cli::parse();
    if let Some(threads) = cli.threads {
        if threads == 0 {
            return Err("--threads must be positive".into());
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|error| error.to_string())?;
    }
    match cli.command {
        Command::Plan => {
            println!(
                "{}",
                serde_json::to_string_pretty(&plan_record()).map_err(|error| error.to_string())?
            );
            Ok(())
        }
        Command::GenerateSegment {
            x_start,
            length,
            batch_steps,
            audit_points,
            audit_seed_x,
            source_commit,
            class_number_receipt_sha256,
            previous_checkpoint,
            output,
            checkpoint_output,
            receipt_output,
        } => generate(
            x_start,
            length,
            batch_steps,
            audit_points,
            audit_seed_x,
            source_commit,
            class_number_receipt_sha256,
            previous_checkpoint.as_deref(),
            &output,
            &checkpoint_output,
            &receipt_output,
        ),
        Command::VerifySegment {
            input,
            checkpoint,
            previous_checkpoint,
            audit_points,
            audit_seed_x,
            receipt_output,
        } => verify(
            &input,
            &checkpoint,
            previous_checkpoint.as_deref(),
            audit_points,
            audit_seed_x,
            &receipt_output,
        ),
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p61-isogeny-class: {error}");
        std::process::exit(1);
    }
}
