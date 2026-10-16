//! Generate and independently replay complete selected-degree p61 adjacency shards.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, ErrorKind, Write};
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::million_store::{digest_file, FileDigest};
use crypto_lib::cryptanalysis::isogeny_walk::p61_adjacency::{
    generate_adjacency_jsonl, plan_record, verify_adjacency_jsonl, AdjacencyConfig,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use serde::{de::DeserializeOwned, Serialize};
use serde_json::{json, Value};

#[derive(Debug, Serialize)]
struct DurableWrite {
    digest: FileDigest,
    sync_status: String,
}

#[derive(Debug, Parser)]
#[command(about = "Generate and replay complete selected-degree p61 adjacency shards")]
struct Cli {
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Print the frozen edge counts, selected degrees, and completion gates.
    Plan,
    /// Expand one immutable window of an admitted vertex segment.
    GenerateShard {
        #[arg(long)]
        source_segment: PathBuf,
        #[arg(long)]
        source_verification_receipt: PathBuf,
        #[arg(long)]
        source_record_offset: u64,
        #[arg(long)]
        source_record_count: u64,
        #[arg(long, default_value_t = 64)]
        batch_vertices: usize,
        #[arg(long)]
        producer_source_commit: String,
        #[arg(long)]
        class_number_receipt_sha256: String,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        receipt_output: PathBuf,
    },
    /// Recompute every root, kernel, codomain, feature, and record byte.
    VerifyShard {
        #[arg(long)]
        source_segment: PathBuf,
        #[arg(long)]
        source_verification_receipt: PathBuf,
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        verifier_source_commit: String,
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

fn verified_source_receipt(
    path: &Path,
    source_segment: &Path,
    source_digest: &FileDigest,
) -> Result<FileDigest, String> {
    let value: Value = read_json(path)?;
    if value.get("status").and_then(Value::as_str) != Some("pass") {
        return Err("source verification receipt is not pass".into());
    }
    let recorded_sha = value
        .pointer("/segment_digest/sha256")
        .and_then(Value::as_str)
        .ok_or("source verification receipt lacks segment_digest.sha256")?;
    if recorded_sha != source_digest.sha256 {
        return Err("source verification receipt binds a different segment digest".into());
    }
    if let Some(recorded_path) = value.get("segment").and_then(Value::as_str) {
        if Path::new(recorded_path) != source_segment {
            return Err("source verification receipt names a different segment path".into());
        }
    }
    digest_file(path)
}

#[allow(clippy::too_many_arguments)]
fn generate(
    source_segment: &Path,
    source_verification_receipt: &Path,
    source_record_offset: u64,
    source_record_count: u64,
    batch_vertices: usize,
    producer_source_commit: String,
    class_number_receipt_sha256: String,
    output: &Path,
    receipt_output: &Path,
) -> Result<(), String> {
    if output.exists() || receipt_output.exists() {
        return Err("refusing to overwrite an adjacency shard or receipt".into());
    }
    let source_digest = digest_file(source_segment)?;
    let source_verification_digest =
        verified_source_receipt(source_verification_receipt, source_segment, &source_digest)?;
    let config = AdjacencyConfig {
        source_record_offset,
        source_record_count,
        batch_vertices,
        producer_source_commit,
        class_number_receipt_sha256,
        source_segment_sha256: source_digest.sha256.clone(),
        source_verification_receipt_sha256: source_verification_digest.sha256.clone(),
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
    let output_file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&partial)
        .map_err(|error| format!("create {}: {error}", partial.display()))?;
    let source_file = File::open(source_segment)
        .map_err(|error| format!("open {}: {error}", source_segment.display()))?;
    let source_reader = BufReader::new(GzDecoder::new(source_file));
    let buffered = BufWriter::new(output_file);
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .operating_system(255)
        .write(buffered, Compression::new(6));
    let result = generate_adjacency_jsonl(source_reader, &mut encoder, &config);
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", partial.display()))?;
    buffered.flush().map_err(|error| error.to_string())?;
    let shard_sync_status = sync_file(buffered.get_ref(), &partial)?;
    let receipt = match result {
        Ok(value) => value,
        Err(error) => {
            return Err(format!(
                "{error}; preserved incomplete adjacency shard at {}",
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
    let adjacency_digest = digest_file(output)?;
    let durable_receipt = json!({
        "status": "pass",
        "source_segment": source_segment,
        "source_segment_digest": source_digest,
        "source_verification_receipt": source_verification_receipt,
        "source_verification_receipt_digest": source_verification_digest,
        "adjacency_shard": output,
        "adjacency_shard_digest": adjacency_digest,
        "durability": { "shard_sync_status": shard_sync_status },
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
    source_segment: &Path,
    source_verification_receipt: &Path,
    input: &Path,
    verifier_source_commit: &str,
    receipt_output: &Path,
) -> Result<(), String> {
    if receipt_output.exists() {
        return Err(format!(
            "refusing to overwrite {}",
            receipt_output.display()
        ));
    }
    let source_digest = digest_file(source_segment)?;
    let source_verification_digest =
        verified_source_receipt(source_verification_receipt, source_segment, &source_digest)?;
    let adjacency_digest = digest_file(input)?;
    let source_file = File::open(source_segment)
        .map_err(|error| format!("open {}: {error}", source_segment.display()))?;
    let adjacency_file =
        File::open(input).map_err(|error| format!("open {}: {error}", input.display()))?;
    let source_reader = BufReader::new(GzDecoder::new(source_file));
    let adjacency_reader = BufReader::new(GzDecoder::new(adjacency_file));
    let receipt = verify_adjacency_jsonl(
        source_reader,
        adjacency_reader,
        verifier_source_commit,
        &source_digest.sha256,
        &source_verification_digest.sha256,
    )?;
    let durable_receipt = json!({
        "status": "pass",
        "source_segment": source_segment,
        "source_segment_digest": source_digest,
        "source_verification_receipt": source_verification_receipt,
        "source_verification_receipt_digest": source_verification_digest,
        "adjacency_shard": input,
        "adjacency_shard_digest": adjacency_digest,
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
        Command::GenerateShard {
            source_segment,
            source_verification_receipt,
            source_record_offset,
            source_record_count,
            batch_vertices,
            producer_source_commit,
            class_number_receipt_sha256,
            output,
            receipt_output,
        } => generate(
            &source_segment,
            &source_verification_receipt,
            source_record_offset,
            source_record_count,
            batch_vertices,
            producer_source_commit,
            class_number_receipt_sha256,
            &output,
            &receipt_output,
        ),
        Command::VerifyShard {
            source_segment,
            source_verification_receipt,
            input,
            verifier_source_commit,
            receipt_output,
        } => verify(
            &source_segment,
            &source_verification_receipt,
            &input,
            &verifier_source_commit,
            &receipt_output,
        ),
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p61-isogeny-adjacency: {error}");
        std::process::exit(1);
    }
}
