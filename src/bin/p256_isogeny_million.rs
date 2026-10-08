//! Generate or independently replay the compact P-256 isogeny-grid certificate.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::million::{
    audit_j_union_paths, generate_jsonl, generate_strip_jsonl, generate_trait_census_jsonl,
    verify_jsonl, verify_strip_jsonl, verify_trait_census_paths, GridConfig, StripConfig,
    PRODUCTION_SIDE,
};
use crypto_lib::cryptanalysis::isogeny_walk::million_store::{
    self, fetch_object, publish_cache, stage_cache, DEFAULT_CHUNK_ROWS,
};
use crypto_lib::cryptanalysis::isogeny_walk::store::S3Loc;
use crypto_lib::cryptanalysis::pollard_collab::cairn::{
    wall_clock, CairnConfig, CairnTransport, Submitter,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use serde_json::json;

#[derive(Debug, Parser)]
#[command(about = "Generate and replay a streaming million-curve P-256 isogeny grid")]
struct Cli {
    /// Rayon worker threads (default: all available cores).
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Generate, audit, and write a deterministic gzip JSONL certificate.
    Generate {
        /// Grid side; production uses 1,000 for exactly 1,000,000 curves.
        #[arg(long, default_value_t = PRODUCTION_SIDE)]
        side: u32,
        /// Rows computed in parallel before their deterministic ordered write.
        #[arg(long)]
        batch_rows: Option<usize>,
        /// Deterministic points per construction-time order audit.
        #[arg(long, default_value_t = 1)]
        audit_points: usize,
        /// Clean 40-hex source commit recorded in the certificate.
        #[arg(long)]
        source_commit: String,
        #[arg(long)]
        output: PathBuf,
    },
    /// Generate a disjoint global-row continuation of the P-256 grid.
    GenerateStrip {
        /// Complete degree-11 row width; production continuation uses 1,000.
        #[arg(long, default_value_t = PRODUCTION_SIDE)]
        width: u32,
        /// First global degree-13 spine row included in this certificate.
        #[arg(long)]
        y_start: u32,
        /// Number of complete rows included in this certificate.
        #[arg(long)]
        height: u32,
        /// Rows computed in parallel before their deterministic ordered write.
        #[arg(long)]
        batch_rows: Option<usize>,
        /// Deterministic points per construction-time order audit.
        #[arg(long, default_value_t = 1)]
        audit_points: usize,
        /// Clean 40-hex source commit recorded in the certificate.
        #[arg(long)]
        source_commit: String,
        #[arg(long)]
        output: PathBuf,
    },
    /// Replay every identity, detector, grid choice, kernel, and isomorphism.
    Verify {
        #[arg(long)]
        input: PathBuf,
        /// Deterministic points per independent order audit.
        #[arg(long, default_value_t = 2)]
        audit_points: usize,
        /// Independent replay starts its order-audit points at this x.
        #[arg(long, default_value_t = 7)]
        audit_seed_x: u64,
        /// Complete rows held for parallel replay at once.
        #[arg(long)]
        batch_rows: Option<usize>,
    },
    /// Independently replay a global-row strip certificate.
    VerifyStrip {
        #[arg(long)]
        input: PathBuf,
        /// Deterministic points per independent order audit.
        #[arg(long, default_value_t = 2)]
        audit_points: usize,
        /// Independent replay starts its order-audit points at this x.
        #[arg(long, default_value_t = 7)]
        audit_seed_x: u64,
        /// Complete rows held for parallel replay at once.
        #[arg(long)]
        batch_rows: Option<usize>,
    },
    /// Check certificate integrity and exact j uniqueness across all inputs.
    AuditJUnion {
        /// Legacy or strip certificates, in deterministic accounting order.
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
    },
    /// Derive a compact structural-trait sidecar from replayed certificates.
    TraitCensus {
        /// Legacy or strip source certificate; repeat in union order.
        #[arg(long = "source", required = true)]
        sources: Vec<PathBuf>,
        /// Clean 40-hex source commit recorded in the census.
        #[arg(long)]
        source_commit: String,
        #[arg(long)]
        output: PathBuf,
    },
    /// Reopen the sources and byte-compare a complete census reconstruction.
    VerifyTraitCensus {
        #[arg(long)]
        input: PathBuf,
        /// Exact source certificate list used by `trait-census`.
        #[arg(long = "source", required = true)]
        sources: Vec<PathBuf>,
    },
    /// Derive deterministic row slices after generation and replay agree.
    StageCache {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        generation_receipt: PathBuf,
        #[arg(long)]
        verification_receipt: PathBuf,
        /// Safe, stable semantic run id used by S3 and Cairn.
        #[arg(long)]
        run: String,
        #[arg(long, default_value_t = DEFAULT_CHUNK_ROWS)]
        chunk_rows: u32,
        #[arg(long)]
        output: PathBuf,
    },
    /// Upload exact bytes and slices, then create the immutable marker last.
    Publish {
        /// S3 root, for example s3://bucket/p256-isogeny-grid.
        #[arg(long)]
        store: String,
        #[arg(long)]
        cache: PathBuf,
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        generation_receipt: PathBuf,
        #[arg(long)]
        verification_receipt: PathBuf,
        #[arg(long)]
        scratch: PathBuf,
        /// Optional Cairn node; S3 remains the byte authority.
        #[arg(long)]
        cairn: Option<String>,
        /// Cairn objective whose checker accepts the coordination receipt.
        #[arg(long)]
        objective: Option<String>,
        /// Stage-0 Cairn nickname. Prefer --identity for signed submissions.
        #[arg(long)]
        submitter: Option<String>,
        /// Identity file created by `cairn identity --out FILE`.
        #[arg(long)]
        identity: Option<PathBuf>,
        /// Persisted commit-reveal outbox; required with --cairn.
        #[arg(long)]
        cairn_state: Option<PathBuf>,
        #[arg(long, default_value_t = 600)]
        cairn_epoch_secs: u64,
    },
    /// Fetch and hash-check the canonical certificate or one cache slice.
    Fetch {
        #[arg(long)]
        store: String,
        #[arg(long)]
        run: String,
        /// Cache part name; omit for the canonical certificate.
        #[arg(long)]
        part: Option<String>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long)]
        scratch: PathBuf,
    },
}

fn default_batch_rows() -> usize {
    (rayon::current_num_threads() * 2).max(1)
}

fn partial_path(output: &Path) -> PathBuf {
    let mut name: OsString = output.as_os_str().to_owned();
    name.push(".partial");
    PathBuf::from(name)
}

fn generate(
    output: &Path,
    side: u32,
    batch_rows: usize,
    audit_points: usize,
    source_commit: String,
) -> Result<(), String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
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
    let result = generate_jsonl(
        &mut encoder,
        &GridConfig {
            side,
            batch_rows,
            audit_points,
            audit_seed_x: 0,
            source_commit,
        },
    );
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", partial.display()))?;
    buffered
        .flush()
        .map_err(|error| format!("flush {}: {error}", partial.display()))?;
    buffered
        .get_ref()
        .sync_all()
        .map_err(|error| format!("sync {}: {error}", partial.display()))?;
    let receipt = match result {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(format!(
                "{error}; preserved the incomplete attempt at {}",
                partial.display()
            ))
        }
    };
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit certificate {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn verify(
    input: &Path,
    audit_points: usize,
    audit_seed_x: u64,
    batch_rows: usize,
) -> Result<(), String> {
    let file = File::open(input).map_err(|error| format!("open {}: {error}", input.display()))?;
    let receipt = if input.extension().and_then(|value| value.to_str()) == Some("gz") {
        verify_jsonl(
            BufReader::new(GzDecoder::new(file)),
            audit_points,
            audit_seed_x,
            batch_rows,
        )?
    } else {
        verify_jsonl(BufReader::new(file), audit_points, audit_seed_x, batch_rows)?
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn generate_strip(
    output: &Path,
    width: u32,
    y_start: u32,
    height: u32,
    batch_rows: usize,
    audit_points: usize,
    source_commit: String,
) -> Result<(), String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
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
    let result = generate_strip_jsonl(
        &mut encoder,
        &StripConfig {
            width,
            y_start,
            height,
            batch_rows,
            audit_points,
            audit_seed_x: 0,
            source_commit,
        },
    );
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", partial.display()))?;
    buffered
        .flush()
        .map_err(|error| format!("flush {}: {error}", partial.display()))?;
    buffered
        .get_ref()
        .sync_all()
        .map_err(|error| format!("sync {}: {error}", partial.display()))?;
    let receipt = match result {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(format!(
                "{error}; preserved the incomplete strip attempt at {}",
                partial.display()
            ))
        }
    };
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit strip certificate {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn verify_strip(
    input: &Path,
    audit_points: usize,
    audit_seed_x: u64,
    batch_rows: usize,
) -> Result<(), String> {
    let file = File::open(input).map_err(|error| format!("open {}: {error}", input.display()))?;
    let receipt = if input.extension().and_then(|value| value.to_str()) == Some("gz") {
        verify_strip_jsonl(
            BufReader::new(GzDecoder::new(file)),
            audit_points,
            audit_seed_x,
            batch_rows,
        )?
    } else {
        verify_strip_jsonl(BufReader::new(file), audit_points, audit_seed_x, batch_rows)?
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn audit_j_union(inputs: &[PathBuf]) -> Result<(), String> {
    let receipt = audit_j_union_paths(inputs)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&receipt).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn encoded_digest(path: &Path) -> Result<million_store::EncodedDigest, String> {
    if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        million_store::digest_gzip(path)
    } else {
        let digest = million_store::digest_file(path)?;
        Ok(million_store::EncodedDigest {
            encoding: "identity".into(),
            stored_sha256: digest.sha256.clone(),
            stored_bytes: digest.bytes,
            content_sha256: digest.sha256,
            content_bytes: digest.bytes,
        })
    }
}

fn trait_census(output: &Path, sources: &[PathBuf], source_commit: String) -> Result<(), String> {
    if output.exists() {
        return Err(format!("refusing to overwrite {}", output.display()));
    }
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
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
    let result = generate_trait_census_jsonl(&mut encoder, sources, source_commit);
    let mut buffered = encoder
        .finish()
        .map_err(|error| format!("finish {}: {error}", partial.display()))?;
    buffered
        .flush()
        .map_err(|error| format!("flush {}: {error}", partial.display()))?;
    buffered
        .get_ref()
        .sync_all()
        .map_err(|error| format!("sync {}: {error}", partial.display()))?;
    let receipt = match result {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(format!(
                "{error}; preserved the incomplete trait census at {}",
                partial.display()
            ))
        }
    };
    fs::rename(&partial, output).map_err(|error| {
        format!(
            "commit trait census {} -> {}: {error}",
            partial.display(),
            output.display()
        )
    })?;
    let artifact = encoded_digest(output)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "receipt": receipt,
            "artifact": {
                "path": output,
                "encoding": artifact.encoding,
                "stored_sha256": artifact.stored_sha256,
                "stored_bytes": artifact.stored_bytes,
                "content_sha256": artifact.content_sha256,
                "content_bytes": artifact.content_bytes,
            }
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn verify_trait_census(input: &Path, sources: &[PathBuf]) -> Result<(), String> {
    let receipt = verify_trait_census_paths(input, sources)?;
    let artifact = encoded_digest(input)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "receipt": receipt,
            "artifact": {
                "path": input,
                "encoding": artifact.encoding,
                "stored_sha256": artifact.stored_sha256,
                "stored_bytes": artifact.stored_bytes,
                "content_sha256": artifact.content_sha256,
                "content_bytes": artifact.content_bytes,
            }
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn stage(
    input: &Path,
    generation_receipt: &Path,
    verification_receipt: &Path,
    run: &str,
    output: &Path,
    chunk_rows: u32,
) -> Result<(), String> {
    let manifest = stage_cache(
        input,
        generation_receipt,
        verification_receipt,
        run,
        output,
        chunk_rows,
    )?;
    let manifest_path = output.join(million_store::CACHE_MANIFEST_FILE);
    let digest = million_store::digest_file(&manifest_path)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "run_id": manifest.run_id,
            "unique_curves": manifest.unique_curves,
            "parts": manifest.parts.len(),
            "canonical": manifest.canonical,
            "manifest_sha256": digest.sha256,
            "manifest_bytes": digest.bytes,
            "output": output,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn publish(
    store_uri: &str,
    cache: &Path,
    input: &Path,
    generation_receipt: &Path,
    verification_receipt: &Path,
    scratch: &Path,
    cairn: Option<String>,
    objective: Option<String>,
    submitter: Option<String>,
    identity: Option<PathBuf>,
    cairn_state: Option<PathBuf>,
    cairn_epoch_secs: u64,
) -> Result<(), String> {
    let loc = S3Loc::parse(store_uri)?;
    let published = publish_cache(
        &loc,
        cache,
        input,
        generation_receipt,
        verification_receipt,
        scratch,
    )?;
    let cairn_status = if let Some(url) = cairn {
        let objective_id = objective.ok_or("--cairn requires --objective")?;
        let state_path = cairn_state.ok_or("--cairn requires --cairn-state")?;
        let who = match (identity, submitter) {
            (Some(path), None) => Submitter::from_identity_file(&path)?,
            (None, Some(name)) => Submitter::Nickname(name),
            (Some(_), Some(_)) => {
                return Err("use either --identity or --submitter, not both".into())
            }
            (None, None) => return Err("--cairn requires --identity or --submitter".into()),
        };
        let mut transport = CairnTransport::open(CairnConfig {
            url,
            objective_id,
            submitter: who,
            answer_objective: None,
            epoch_secs: cairn_epoch_secs,
            state_path: Some(state_path),
            clock: wall_clock(),
        })?;
        let (revealed, refused) = transport.reveal_pending()?;
        let report = transport.publish_artifact(
            published.cairn_receipt.artifact()?,
            published.cairn_receipt.claim_key(),
        )?;
        json!({
            "committed": report.committed,
            "skipped": report.skipped,
            "revealed": revealed,
            "refused": refused,
            "pending": transport.pending(),
            "submitter": transport.submitter(),
        })
    } else {
        if objective.is_some() || submitter.is_some() || identity.is_some() || cairn_state.is_some()
        {
            return Err("Cairn options require --cairn".into());
        }
        serde_json::Value::Null
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "marker_uri": published.marker_uri,
            "marker_sha256": published.marker_sha256,
            "marker_created": published.created,
            "canonical": published.marker.canonical,
            "parts": published.marker.parts.len(),
            "cairn_receipt": published.cairn_receipt,
            "cairn": cairn_status,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn fetch(
    store_uri: &str,
    run: &str,
    part: Option<&str>,
    output: &Path,
    scratch: &Path,
) -> Result<(), String> {
    let loc = S3Loc::parse(store_uri)?;
    let object = fetch_object(&loc, run, part, output, scratch)?;
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "status": "pass",
            "run_id": run,
            "part": part,
            "output": output,
            "object": object,
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
        Command::Generate {
            side,
            batch_rows,
            audit_points,
            source_commit,
            output,
        } => generate(
            &output,
            side,
            batch_rows.unwrap_or_else(default_batch_rows),
            audit_points,
            source_commit,
        ),
        Command::GenerateStrip {
            width,
            y_start,
            height,
            batch_rows,
            audit_points,
            source_commit,
            output,
        } => generate_strip(
            &output,
            width,
            y_start,
            height,
            batch_rows.unwrap_or_else(default_batch_rows),
            audit_points,
            source_commit,
        ),
        Command::Verify {
            input,
            audit_points,
            audit_seed_x,
            batch_rows,
        } => verify(
            &input,
            audit_points,
            audit_seed_x,
            batch_rows.unwrap_or_else(default_batch_rows),
        ),
        Command::VerifyStrip {
            input,
            audit_points,
            audit_seed_x,
            batch_rows,
        } => verify_strip(
            &input,
            audit_points,
            audit_seed_x,
            batch_rows.unwrap_or_else(default_batch_rows),
        ),
        Command::AuditJUnion { inputs } => audit_j_union(&inputs),
        Command::TraitCensus {
            sources,
            source_commit,
            output,
        } => trait_census(&output, &sources, source_commit),
        Command::VerifyTraitCensus { input, sources } => verify_trait_census(&input, &sources),
        Command::StageCache {
            input,
            generation_receipt,
            verification_receipt,
            run,
            chunk_rows,
            output,
        } => stage(
            &input,
            &generation_receipt,
            &verification_receipt,
            &run,
            &output,
            chunk_rows,
        ),
        Command::Publish {
            store,
            cache,
            input,
            generation_receipt,
            verification_receipt,
            scratch,
            cairn,
            objective,
            submitter,
            identity,
            cairn_state,
            cairn_epoch_secs,
        } => publish(
            &store,
            &cache,
            &input,
            &generation_receipt,
            &verification_receipt,
            &scratch,
            cairn,
            objective,
            submitter,
            identity,
            cairn_state,
            cairn_epoch_secs,
        ),
        Command::Fetch {
            store,
            run,
            part,
            output,
            scratch,
        } => fetch(&store, &run, part.as_deref(), &output, &scratch),
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p256-isogeny-million: {error}");
        std::process::exit(1);
    }
}
