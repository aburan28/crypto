//! Generate or independently replay the compact P-256 isogeny-grid certificate.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::million::{
    generate_jsonl, verify_jsonl, GridConfig, PRODUCTION_SIDE,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};

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
    }
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p256-isogeny-million: {error}");
        std::process::exit(1);
    }
}
