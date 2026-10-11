//! Frozen exact CM class-relation experiment for NIST P-192.

use crypto_lib::cryptanalysis::p192_cm_relation_search as experiment;

use clap::{Parser, Subcommand};
use experiment::{
    build_generator_report, build_order_report, build_pocklington_certificate, run_exact_search,
    unsupported_weights, verify_generator_report, verify_order_report, verify_search_report,
    verify_unsupported_weights, GeneratorReport, OrderReport, SearchReport, WeightsReport,
    HARD_STATE_CAP,
};
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(about = "Exact weighted CM class-relation census for the valid NIST P-192 isogeny class")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Certify the pinned curve/order facts and prove the large D cofactor prime.
    ProveOrder {
        #[arg(long)]
        curve: String,
        #[arg(long)]
        pocklington: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Recompute the thirteen frozen oriented prime-ideal generators.
    Generators {
        #[arg(long)]
        order: PathBuf,
        #[arg(long)]
        max_ell: u64,
        #[arg(long)]
        include: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Explicitly record the currently unavailable HOT/COLD map calibration.
    Calibrate {
        #[arg(long)]
        generators: PathBuf,
        #[arg(long)]
        endpoints: u64,
        #[arg(long)]
        warmup: u64,
        #[arg(long)]
        reps: u64,
        #[arg(long)]
        seed: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Exhaust the frozen exact norm box.  This is the long (about 5M-state) lane.
    Search {
        #[arg(long)]
        mode: String,
        #[arg(long)]
        max_ell: u64,
        #[arg(long)]
        max_degree_bits: u32,
        #[arg(long)]
        max_states: u64,
        #[arg(long)]
        canonical_conjugates: bool,
        #[arg(long)]
        generators: PathBuf,
        #[arg(long)]
        weights: PathBuf,
        #[arg(long)]
        word_radius_control: u32,
        #[arg(long)]
        out: PathBuf,
    },
    /// Independently replay a completed search certificate.
    Verify {
        #[arg(long)]
        certificate: PathBuf,
        #[arg(long)]
        independent_ideal_hnf: bool,
        #[arg(long)]
        replay_maps: bool,
        #[arg(long)]
        out: PathBuf,
    },
}

fn read_json<T: DeserializeOwned>(path: &Path) -> Result<T, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes).map_err(|error| format!("parse {}: {error}", path.display()))
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<(), String> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    bytes.push(b'\n');
    fs::write(path, bytes).map_err(|error| format!("write {}: {error}", path.display()))
}

fn require_p192(curve: &str) -> Result<(), String> {
    if curve.eq_ignore_ascii_case("p192")
        || curve.eq_ignore_ascii_case("p-192")
        || curve.eq_ignore_ascii_case("secp192r1")
    {
        Ok(())
    } else {
        Err(format!(
            "unsupported curve `{curve}`; this frozen binary accepts only P-192"
        ))
    }
}

fn run(command: Command) -> Result<(), String> {
    match command {
        Command::ProveOrder {
            curve,
            pocklington,
            out,
        } => {
            require_p192(&curve)?;
            let certificate = build_pocklington_certificate()?;
            let report = build_order_report(&certificate)?;
            verify_order_report(&report, &certificate)?;
            write_json(&pocklington, &certificate)?;
            write_json(&out, &report)?;
            eprintln!("P-192 order/CM certificate: {}", report.status);
        }
        Command::Generators {
            order,
            max_ell,
            include,
            out,
        } => {
            let normalized: std::collections::BTreeSet<_> =
                include.split(',').map(str::trim).collect();
            if normalized != std::collections::BTreeSet::from(["ramified", "split"]) {
                return Err("the frozen generator set requires --include split,ramified".to_owned());
            }
            let order: OrderReport = read_json(&order)?;
            let pocklington = build_pocklington_certificate()?;
            verify_order_report(&order, &pocklington)?;
            let report = build_generator_report(&order, max_ell)?;
            verify_generator_report(&report)?;
            write_json(&out, &report)?;
            eprintln!(
                "P-192 generators: {} records ({})",
                report.generators.len(),
                report.status
            );
        }
        Command::Calibrate {
            generators,
            endpoints,
            warmup,
            reps,
            seed,
            out,
        } => {
            let generators: GeneratorReport = read_json(&generators)?;
            verify_generator_report(&generators)?;
            if endpoints != 16 || warmup != 1 || reps != 11 || seed != "0x50313932434d5231" {
                return Err("calibration arguments differ from the frozen protocol".to_owned());
            }
            let report = unsupported_weights(&seed, endpoints, reps, warmup);
            write_json(&out, &report)?;
            eprintln!("P-192 HOT/COLD calibration: {}", report.status);
        }
        Command::Search {
            mode,
            max_ell,
            max_degree_bits,
            max_states,
            canonical_conjugates,
            generators,
            weights,
            word_radius_control,
            out,
        } => {
            if mode != "exact-norm"
                || max_ell != 113
                || max_degree_bits != 48
                || max_states != HARD_STATE_CAP
                || !canonical_conjugates
            {
                return Err(
                    "search arguments differ from the frozen exact-norm protocol".to_owned(),
                );
            }
            let generators: GeneratorReport = read_json(&generators)?;
            let weights: WeightsReport = read_json(&weights)?;
            verify_unsupported_weights(&weights)?;
            let report =
                run_exact_search(&generators, Some(&weights), max_states, word_radius_control)?;
            write_json(&out, &report)?;
            eprintln!(
                "P-192 exact census: {} canonical vectors, {} non-scalar hits",
                report.enumeration_forward.canonical_vectors_including_zero,
                report.non_scalar_relation_count
            );
        }
        Command::Verify {
            certificate,
            independent_ideal_hnf,
            replay_maps,
            out,
        } => {
            let report: SearchReport = read_json(&certificate)?;
            let verification = verify_search_report(&report, independent_ideal_hnf, replay_maps)?;
            write_json(&out, &verification)?;
            eprintln!("P-192 relation verification: {}", verification.verdict);
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse().command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_alias_gate_is_narrow() {
        assert!(require_p192("p192").is_ok());
        assert!(require_p192("secp192r1").is_ok());
        assert!(require_p192("p256").is_err());
    }

    #[test]
    fn report_types_round_trip() {
        let certificate: experiment::PocklingtonCertificate =
            build_pocklington_certificate().unwrap();
        let text = serde_json::to_string(&certificate).unwrap();
        let decoded: experiment::PocklingtonCertificate = serde_json::from_str(&text).unwrap();
        experiment::verify_pocklington_certificate(&decoded).unwrap();
    }
}
