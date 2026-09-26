//! `ecc2k-guard` — checks the ECC2K-130 rho campaign can be trusted with
//! its budget.
//!
//! ```text
//! ecc2k-guard cycles                              # no fruitless cycle up to 8 steps
//! ecc2k-guard cycles --max-length 16 --json cert.json
//! ecc2k-guard seed-generate --out campaign-seed.json
//! ecc2k-guard seed-verify campaign-seed.json
//! ecc2k-guard seed-decode --seed 00c4005dffff0007
//! ```
//!
//! Every subcommand exits non-zero when its check fails.

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::ecc2k130_guard::{
    self as guard, ChallengeAnchor, CycleVerdict, SeedProvenance,
};

#[derive(Parser)]
#[command(name = "ecc2k-guard", about = "ECC2K-130 rho campaign guards")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Rule out every fruitless cycle of the shipping walk up to a length.
    Cycles {
        /// Longest cycle to rule out, in steps.
        #[arg(long, default_value_t = 8)]
        max_length: usize,
        /// Field degrees to check; the client's five by default.
        #[arg(long, value_delimiter = ',')]
        curves: Vec<u32>,
        /// Also write the certificate here.
        #[arg(long)]
        json: Option<PathBuf>,
    },
    /// Draw a campaign master seed and write its provenance record.
    SeedGenerate {
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value = "ecc2k-130")]
        campaign: String,
        /// Use this recorded public value (at least 32 bytes, hex) instead
        /// of the OS CSPRNG, e.g. a block hash published after the
        /// commitment was announced.
        #[arg(long, requires = "beacon_label")]
        beacon_hex: Option<String>,
        /// Where the beacon value came from, stored verbatim.
        #[arg(long)]
        beacon_label: Option<String>,
    },
    /// Re-derive a provenance record and check every field.
    SeedVerify { path: PathBuf },
    /// Show what the live client derives from a walk seed.
    SeedDecode {
        /// A seed as recorded in a distinguished-point record, hex.
        #[arg(long, conflicts_with_all = ["run_id", "walk"])]
        seed: Option<String>,
        #[arg(long, requires = "walk")]
        run_id: Option<u16>,
        #[arg(long, requires = "run_id")]
        walk: Option<u32>,
        #[arg(long, default_value_t = 0)]
        restart: u16,
    },
}

fn cycles(max_length: usize, curves: Vec<u32>, json: Option<PathBuf>) -> Result<bool, String> {
    let curves = if curves.is_empty() {
        guard::CURVES.to_vec()
    } else {
        curves
    };
    let cert = guard::certify(&curves, max_length);
    println!(
        "no-fruitless-cycle check, steps j in {:?}, cycles up to length {max_length}",
        cert.step_exponents
    );
    for c in &cert.curves {
        let expected = c.expected_by_chance;
        match c.verdict {
            CycleVerdict::Skipped => println!(
                "  m={:<4} skipped: {}",
                c.m,
                c.skip_reason.as_deref().unwrap_or("")
            ),
            CycleVerdict::Clean => println!(
                "  m={:<4} clean: {} multisets, none hit an orbit ({expected:.2e} expected by chance)",
                c.m, c.multisets_checked
            ),
            CycleVerdict::Coincidence | CycleVerdict::Real => {
                println!(
                    "  m={:<4} {} hit(s) in {} multisets, {} ({expected:.2e} expected by chance)",
                    c.m,
                    c.hits.len(),
                    c.multisets_checked,
                    if c.verdict == CycleVerdict::Real {
                        "REAL"
                    } else {
                        "coincidence at this size"
                    },
                );
                for h in c.hits.iter().take(2) {
                    println!(
                        "         length {}: j = {:?} closes on {}",
                        h.length, h.exponents, h.orbit_scalar
                    );
                }
            }
            CycleVerdict::Degenerate => {
                for why in &c.degenerate {
                    println!("  m={:<4} DEGENERATE: {why}", c.m);
                }
            }
        }
        match c.matches_challenge_constants {
            Some(true) => println!(
                "  m={:<4} ell and lambda equal the generated challenge constants",
                c.m
            ),
            Some(false) => println!(
                "  m={:<4} MISMATCH: ell or lambda differs from the generated challenge constants",
                c.m
            ),
            None => {}
        }
    }
    if let Some(path) = json {
        let body = serde_json::to_string_pretty(&cert).map_err(|e| e.to_string())? + "\n";
        std::fs::write(&path, body).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    if cert.passed {
        println!("\nNo cycle mechanism: the iteration is equivariant, and no step multiset up to the bound closes an orbit.");
    } else {
        println!("\nA real hit or a degenerate step means a walk can loop without reaching a distinguished point: the step schedule has to change.");
    }
    Ok(cert.passed)
}

fn seed_generate(
    out: PathBuf,
    campaign: String,
    beacon_hex: Option<String>,
    beacon_label: Option<String>,
) -> Result<bool, String> {
    let (entropy, source) = match beacon_hex {
        Some(h) => {
            let bytes = hex::decode(h.trim()).map_err(|e| format!("--beacon-hex: {e}"))?;
            if bytes.len() < 32 {
                return Err(format!(
                    "--beacon-hex carries {} bytes; at least 32 are required",
                    bytes.len()
                ));
            }
            let label = beacon_label.unwrap_or_default();
            (bytes, format!("beacon: {label}"))
        }
        None => (
            guard::os_entropy()?.to_vec(),
            "getrandom: OS CSPRNG".to_string(),
        ),
    };
    let created = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_secs();
    let record = guard::generate_seed(
        ChallengeAnchor::ecc2k130(&campaign),
        &entropy,
        &source,
        created,
    );
    guard::verify_seed(&record).map_err(|f| f.join("; "))?;
    if out.exists() {
        return Err(format!(
            "{} exists; a campaign seed is written once and never overwritten",
            out.display()
        ));
    }
    let body = serde_json::to_string_pretty(&record).map_err(|e| e.to_string())? + "\n";
    std::fs::write(&out, body).map_err(|e| format!("{}: {e}", out.display()))?;
    println!("wrote {}", out.display());
    println!("  entropy        {source}");
    println!("  master seed    {}", record.master_seed_hex);
    println!("  commitment     {}", record.commitment_blake3);
    Ok(true)
}

fn seed_verify(path: PathBuf) -> Result<bool, String> {
    let bytes = std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))?;
    let record: SeedProvenance =
        serde_json::from_slice(&bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    match guard::verify_seed(&record) {
        Ok(()) => {
            println!(
                "{}: verified; master seed {} re-derives from the challenge anchor and the recorded entropy ({})",
                path.display(),
                record.master_seed_hex,
                record.entropy_source
            );
            Ok(true)
        }
        Err(failures) => {
            for f in failures {
                println!("{}: FAILED: {f}", path.display());
            }
            Ok(false)
        }
    }
}

fn seed_decode(
    seed: Option<String>,
    run_id: Option<u16>,
    walk: Option<u32>,
    restart: u16,
) -> Result<bool, String> {
    let seed = match (seed, run_id, walk) {
        (Some(s), _, _) => u64::from_str_radix(s.trim().trim_start_matches("0x"), 16)
            .map_err(|e| format!("--seed: {e}"))?,
        (None, Some(r), Some(w)) => guard::legacy_seed(r, w, restart),
        _ => return Err("give --seed HEX, or --run-id and --walk".into()),
    };
    let (run, walk, restarts) = guard::decode_legacy_seed(seed);
    let [c0, c1] = guard::legacy_start_bits(seed);
    println!("seed          {seed:016x}");
    println!("run id        {run}");
    println!("walk index    {walk}");
    println!("restarts      {restarts} of 65535 before the client refuses to wrap");
    println!(
        "start bits    c0 {c0:016x}  c1 {c1:016x}  (splitmix64; start = Q + sum c_i sigma^i(P))"
    );
    Ok(true)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let outcome = match cli.command {
        Command::Cycles {
            max_length,
            curves,
            json,
        } => cycles(max_length, curves, json),
        Command::SeedGenerate {
            out,
            campaign,
            beacon_hex,
            beacon_label,
        } => seed_generate(out, campaign, beacon_hex, beacon_label),
        Command::SeedVerify { path } => seed_verify(path),
        Command::SeedDecode {
            seed,
            run_id,
            walk,
            restart,
        } => seed_decode(seed, run_id, walk, restart),
    };
    match outcome {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::from(2)
        }
    }
}
