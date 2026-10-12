//! `ecc2k-dp-ingest` — CLI stub for the ECC2K-130 distinguished-point ingest.
//!
//! Phase A of the Rust port: the pure library is in place and tested. This
//! binary accepts the same flags as `ecc2k130/aws/dp_ingest.py` so deploy
//! wrappers can start pointing at it, but every mode that needs Postgres or
//! S3 exits with a clear message. The live daemon remains Python until the
//! parity harness of Phase B lands. See
//! `research/ecc2k130_dp_ingest_rust_20261007/PROTOCOL.md`.

use clap::Parser;
use ecc2k130_dp_ingest::{decode, is_witness_v2, table_record_body, DP_MAGIC_V2, RECORD_BYTES};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser, Debug)]
#[command(
    name = "ecc2k-dp-ingest",
    about = "ECC2K-130 distinguished-point ingest (pure core; Python still runs live)"
)]
struct Cli {
    /// Campaign S3 bucket holding dp/
    #[arg(long, env = "RHO_BUCKET")]
    bucket: Option<String>,
    /// Bucket that receives status.json
    #[arg(long, env = "RHO_STATUS_BUCKET")]
    status_bucket: Option<String>,
    /// Campaign id
    #[arg(long, env = "RHO_CAMPAIGN", default_value = "ecc2k-130")]
    campaign: String,
    /// Drain one pass and exit
    #[arg(long)]
    once: bool,
    /// Report the backlog and exit
    #[arg(long)]
    pending: bool,
    /// Re-check already-stored points against a corpus
    #[arg(long)]
    verify: bool,
    /// Read checkpoints only; no database
    #[arg(long)]
    work: bool,
    /// Rebuild rollup counters from the points table
    #[arg(long)]
    recount: bool,
    /// Skip CREATE INDEX CONCURRENTLY (Lambda default)
    #[arg(long)]
    no_index: bool,
    /// Ingest threads
    #[arg(long, env = "RHO_INGEST_THREADS", default_value_t = 4)]
    threads: u32,
    /// CloudWatch namespace for ingest metrics
    #[arg(long, env = "RHO_METRIC_NAMESPACE")]
    metric_namespace: Option<String>,
    /// Offline: decode a local .bin and print JSON lines (Phase A)
    #[arg(long)]
    decode_file: Option<PathBuf>,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    if let Some(path) = cli.decode_file {
        return decode_file(&path);
    }
    eprintln!(
        "ecc2k-dp-ingest: Phase A stub — pure decode/cutoff/envelope logic is in the library; \
         Postgres and S3 land in Phase B. The live daemon is still \
         ecc2k130/aws/dp_ingest.py. See research/ecc2k130_dp_ingest_rust_20261007/PROTOCOL.md."
    );
    eprintln!(
        "ecc2k-dp-ingest: refusing to run {} against bucket={:?} (use --decode-file for the offline path)",
        mode(&cli),
        cli.bucket
    );
    let _ = (
        cli.status_bucket,
        cli.campaign,
        cli.threads,
        cli.metric_namespace,
        cli.no_index,
    );
    ExitCode::from(2)
}

fn mode(cli: &Cli) -> &'static str {
    if cli.verify {
        "verify"
    } else if cli.pending {
        "pending"
    } else if cli.work {
        "work"
    } else if cli.recount {
        "recount"
    } else if cli.once {
        "once"
    } else {
        "poll"
    }
}

fn decode_file(path: &PathBuf) -> ExitCode {
    let body = match std::fs::read(path) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("ecc2k-dp-ingest: {}: {e}", path.display());
            return ExitCode::FAILURE;
        }
    };
    if is_witness_v2(&body) {
        eprintln!(
            "ecc2k-dp-ingest: {}: a WITNESS=1 corpus ({}…); this store takes 32-byte records only",
            path.display(),
            String::from_utf8_lossy(DP_MAGIC_V2)
        );
        return ExitCode::from(3);
    }
    let records = match table_record_body(&body) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("ecc2k-dp-ingest: {}: {e}", path.display());
            return ExitCode::FAILURE;
        }
    };
    if !records.len().is_multiple_of(RECORD_BYTES) {
        eprintln!(
            "ecc2k-dp-ingest: {}: length {} is not a multiple of {RECORD_BYTES}",
            path.display(),
            records.len()
        );
        return ExitCode::FAILURE;
    }
    for chunk in records.as_chunks::<RECORD_BYTES>().0 {
        let row = decode(chunk);
        println!(
            "{{\"point_key\":\"{}\",\"a\":\"{}\",\"walk_seed\":\"{}\"}}",
            hex::encode(row.point_key),
            hex::encode(row.a),
            hex::encode(row.walk_seed)
        );
    }
    ExitCode::SUCCESS
}
