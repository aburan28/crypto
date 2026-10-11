//! `ecc2k-status` — the ECC2K-130 campaign dashboard from the slot table:
//! who is walking, how fast, how far, and how long their walks run.
//!
//! ```text
//! ecc2k-status                 one snapshot
//! ecc2k-status --watch 60      refresh every 60 s
//! ecc2k-status --json          machine-readable
//! ```
//!
//! The slots come from `--local DIR` (a rehearsal store), `--table`
//! (DynamoDB) or `--bucket` (the `slots/` objects), preferred in that
//! order, each defaulting to `ECC_LOCAL_STORE`, `ECC_TABLE` and
//! `ECC_BUCKET`. The AWS sources run the `aws` CLI. See
//! `crypto_lib::cryptanalysis::ecc2k130_status` for what is counted.

use std::io::{self, Write};
use std::process::ExitCode;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use clap::Parser;
use crypto_lib::cryptanalysis::ecc2k130_pyjson::float_literal;
use crypto_lib::cryptanalysis::ecc2k130_status::{parse_walks, report, scan, Source, PRESET_WALKS};

#[derive(Parser)]
#[command(
    name = "ecc2k-status",
    about = "ECC2K-130 campaign dashboard from the slot table"
)]
struct Cli {
    /// Campaign bucket (slots/ prefix) [default: $ECC_BUCKET]
    #[arg(long)]
    bucket: Option<String>,
    /// DynamoDB table, if the slots live there [default: $ECC_TABLE]
    #[arg(long)]
    table: Option<String>,
    /// Rehearsal store directory [default: $ECC_LOCAL_STORE]
    #[arg(long)]
    local: Option<String>,
    /// Walks assumed for a slot that does not report its own (workers x batch)
    #[arg(long, default_value_t = PRESET_WALKS, value_parser = parse_walks,
          allow_negative_numbers = true)]
    walks: i128,
    /// Refresh every this many seconds; 0 prints one snapshot
    #[arg(long, default_value_t = 0.0, value_parser = parse_watch)]
    watch: f64,
    /// Print the summary as JSON
    #[arg(long)]
    json: bool,
}

fn parse_watch(s: &str) -> Result<f64, String> {
    match float_literal(s) {
        Some(x) if x.is_finite() && x >= 0.0 => Ok(x),
        _ => Err(format!("invalid watch interval {s:?}: seconds, 0 or more")),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let pick = |flag: Option<String>, var: &str| {
        flag.or_else(|| std::env::var(var).ok())
            .filter(|v| !v.is_empty())
    };
    let source = match (
        pick(cli.local, "ECC_LOCAL_STORE"),
        pick(cli.table, "ECC_TABLE"),
        pick(cli.bucket, "ECC_BUCKET"),
    ) {
        (Some(local), _, _) => Source::Local(local.into()),
        (None, Some(table), _) => Source::Table(table),
        (None, None, Some(bucket)) => Source::Bucket(bucket),
        (None, None, None) => {
            eprintln!("give --bucket (or ECC_BUCKET), --table or --local");
            return ExitCode::FAILURE;
        }
    };
    loop {
        let mut out = String::new();
        let result = scan(&source, "aws").and_then(|items| {
            let now = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_or(0.0, |d| d.as_secs_f64());
            report(&mut out, &items, cli.walks, now, cli.json)
        });
        print!("{out}");
        let _ = io::stdout().flush();
        if let Err(e) = &result {
            eprintln!("status failed: {e}");
        }
        if cli.watch == 0.0 {
            return if result.is_ok() {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            };
        }
        std::thread::sleep(Duration::from_secs_f64(cli.watch));
        println!();
    }
}
