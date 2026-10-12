//! Generate and replay the bounded native P-256 degree-11 walk certificate.

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::p256_isogeny_walk::{
    from_json_lines, generate, to_json_lines, verify, verify_independent_start_fixture,
    DEFAULT_STEPS,
};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(about = "Generate or verify the bounded P-256 degree-11 walk certificate")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Generate the deterministic walk, replay it, then write its certificate.
    Generate {
        #[arg(long, default_value_t = DEFAULT_STEPS)]
        steps: u64,
        #[arg(long)]
        output: PathBuf,
    },
    /// Replay a previously generated certificate.
    Verify {
        #[arg(long)]
        input: PathBuf,
    },
    /// Run the independent division-polynomial/Vélu start-edge check.
    Fixture,
}

fn read_certificate(path: &Path) -> Result<Vec<u8>, String> {
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    let mut bytes = Vec::new();
    if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        GzDecoder::new(file)
            .read_to_end(&mut bytes)
            .map_err(|error| format!("decompress {}: {error}", path.display()))?;
    } else {
        let mut file = file;
        file.read_to_end(&mut bytes)
            .map_err(|error| format!("read {}: {error}", path.display()))?;
    }
    Ok(bytes)
}

fn write_certificate(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let file = File::create(path).map_err(|error| format!("create {}: {error}", path.display()))?;
    if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        let mut encoder = GzEncoder::new(file, Compression::best());
        encoder
            .write_all(bytes)
            .map_err(|error| format!("compress {}: {error}", path.display()))?;
        encoder
            .finish()
            .map_err(|error| format!("finish {}: {error}", path.display()))?;
    } else {
        let mut file = file;
        file.write_all(bytes)
            .map_err(|error| format!("write {}: {error}", path.display()))?;
    }
    Ok(())
}

fn run() -> Result<(), String> {
    match Cli::parse().command {
        Command::Generate { steps, output } => {
            let certificate = generate(steps)?;
            let bytes = to_json_lines(&certificate)?;
            write_certificate(&output, &bytes)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&certificate.summary)
                    .map_err(|error| error.to_string())?
            );
        }
        Command::Verify { input } => {
            let bytes = read_certificate(&input)?;
            let mut certificate = from_json_lines(&bytes)?;
            verify(&mut certificate)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&certificate.summary)
                    .map_err(|error| error.to_string())?
            );
        }
        Command::Fixture => {
            let fixtures = verify_independent_start_fixture()?;
            println!(
                "{}",
                serde_json::to_string_pretty(&fixtures).map_err(|error| error.to_string())?
            );
        }
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("p256-isogeny-walk: {error}");
        std::process::exit(1);
    }
}
