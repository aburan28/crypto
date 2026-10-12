//! `ecc2k-merge` — the ECC2K-130 distinguished-point merge, byte-for-byte
//! with `ecc2k130/aws/merge.py` and taking the same arguments.
//!
//! ```text
//! ecc2k-merge --legacy --work WORK --local DIR --curve 41 --dp-weight 18 --client ecc2k130/ecc2k130-cpu
//! ecc2k-merge --campaign campaign.json --work /data/merge-v2 --s3 s3://ecc2k130-<account>/dp/ --client ecc2k130/ecc2k130-cpu
//! ```
//!
//! Exit status: 0 done (or solved), 2 collisions still unsolved or a usage
//! error, 1 a failure that leaves the merge state as it was.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{ArgGroup, Parser};
#[cfg(unix)]
use crypto_lib::cryptanalysis::ecc2k130_merge::{run, MergeArgs, MergeError};

#[derive(Parser)]
#[command(
    name = "ecc2k-merge",
    about = "Merge distinguished-point corpora from many workers and solve any collision"
)]
#[command(group(ArgGroup::new("source").required(true).args(["local", "s3"])))]
struct Cli {
    /// State directory (buckets, offsets, results).
    #[arg(long)]
    work: PathBuf,
    /// Directory of *.bin corpus files.
    #[arg(long)]
    local: Option<PathBuf>,
    /// s3://bucket/prefix/ of corpus deltas to sync.
    #[arg(long)]
    s3: Option<String>,
    /// Host client used to rewalk and solve.
    #[arg(long, default_value = concat!(env!("CARGO_MANIFEST_DIR"), "/ecc2k130/ecc2k130-cpu"))]
    client: PathBuf,
    #[arg(long, default_value_t = 131, allow_negative_numbers = true)]
    curve: i128,
    #[arg(long, default_value = "sigma", value_parser = ["sigma", "table"])]
    walk: String,
    /// Bucket count for a new work dir.
    #[arg(long, default_value_t = 4096, allow_negative_numbers = true)]
    buckets: i128,
    #[arg(long)]
    detect_only: bool,
    /// Strict campaign.json; require committed checksummed chunks.
    #[arg(long)]
    campaign: Option<PathBuf>,
    /// Explicitly allow unversioned raw corpora (not certified).
    #[arg(long)]
    legacy: bool,
    #[arg(long, default_value_t = 3600.0, allow_negative_numbers = true)]
    solve_timeout: f64,
    /// The campaign's cutoff; the rewalk must stop at the same points the
    /// workers reported.
    #[arg(long, allow_negative_numbers = true)]
    dp_weight: Option<i128>,
    /// Extra client argument for the solve step (repeatable).
    #[arg(long, allow_hyphen_values = true)]
    client_arg: Vec<String>,
}

#[cfg(not(unix))]
fn main() -> ExitCode {
    eprintln!("ecc2k-merge requires a unix host");
    ExitCode::from(1)
}

#[cfg(unix)]
fn main() -> ExitCode {
    let cli = Cli::parse();
    let args = MergeArgs {
        work: cli.work,
        local: cli.local,
        s3: cli.s3,
        client: cli.client,
        curve: cli.curve,
        walk: cli.walk,
        buckets: cli.buckets,
        detect_only: cli.detect_only,
        campaign: cli.campaign,
        legacy: cli.legacy,
        solve_timeout: cli.solve_timeout,
        dp_weight: cli.dp_weight,
        client_arg: cli.client_arg,
    };
    let mut stdout = std::io::stdout().lock();
    match run(args, &mut stdout) {
        Ok(code) => ExitCode::from(code as u8),
        Err(e) => {
            match &e {
                MergeError::Usage(m) => eprintln!("ecc2k-merge: error: {m}"),
                MergeError::Fatal(m) => eprintln!("ecc2k-merge: {m}"),
            }
            ExitCode::from(e.exit_code() as u8)
        }
    }
}
