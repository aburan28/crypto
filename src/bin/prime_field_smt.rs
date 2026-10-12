use clap::{Parser, Subcommand, ValueEnum};
use crypto_lib::cryptanalysis::prime_field_smt::{
    blake3_file, emit_finite_field_s3_smt2, emit_finite_field_smt2, emit_smt2,
    exhaustive_reference, lift_s3_abscissae, parse_abscissa_output_for_instance,
    parse_solver_output_for_instance, run_cvc5_with_options, validate_instance, verify_witness,
    InstanceValidation, ParsedAbscissaOutput, ParsedSolverOutput, PrimeFieldSmtInstance,
    VerificationReceipt,
};
use serde::Serialize;
use serde_json::json;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Duration;

#[derive(Parser)]
#[command(
    name = "prime-field-smt",
    about = "Exact SMT two-point decomposition over a prime field"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum EncodingKind {
    /// Portable bit-vectors with widened arithmetic and explicit remainder.
    BitVector,
    /// cvc5 native prime finite fields with Boolean field-bit bounds.
    FiniteField,
    /// Native finite-field Semaev S3 with exact liftable-x membership.
    FiniteFieldS3,
}

#[derive(Subcommand)]
enum Command {
    /// Validate an instance without emitting or solving it.
    Check {
        #[arg(long)]
        instance: PathBuf,
    },
    /// Emit SMT-LIB and a content-addressed encoding receipt.
    Emit {
        #[arg(long)]
        instance: PathBuf,
        #[arg(long)]
        smt2: PathBuf,
        #[arg(long)]
        receipt: PathBuf,
        #[arg(long, value_enum, default_value_t = EncodingKind::BitVector)]
        encoding: EncodingKind,
    },
    /// Solve exactly by enumerating the factor set and testing R-P membership.
    Reference {
        #[arg(long)]
        instance: PathBuf,
        #[arg(long, default_value_t = 1_000_000)]
        max_factor_points: usize,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Parse and independently verify an existing solver stdout file.
    Verify {
        #[arg(long)]
        instance: PathBuf,
        #[arg(long)]
        solver_output: PathBuf,
        #[arg(long)]
        out: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = EncodingKind::BitVector)]
        encoding: EncodingKind,
    },
    /// Emit, run cvc5, retain raw output, and verify any SAT model.
    Solve {
        #[arg(long)]
        instance: PathBuf,
        #[arg(long)]
        solver: PathBuf,
        #[arg(long, default_value_t = 120_000)]
        timeout_ms: u64,
        #[arg(long)]
        out_dir: PathBuf,
        #[arg(long, value_enum, default_value_t = EncodingKind::BitVector)]
        encoding: EncodingKind,
    },
}

#[derive(Serialize)]
struct ByteReceipt {
    bytes: u64,
    blake3: String,
}

#[derive(Serialize)]
struct SolverReceipt {
    executable: String,
    blake3_before: String,
    blake3_after: String,
    unchanged_during_run: bool,
    arguments: Vec<String>,
    requested_timeout_ms: u64,
    harness_grace_ms: u64,
    elapsed_ms: u64,
    timed_out: bool,
    success: bool,
    exit_code: Option<i32>,
    stdout: ByteReceipt,
    stderr: ByteReceipt,
}

#[derive(Serialize)]
struct SolveReceipt {
    schema: String,
    input_path: String,
    input: ByteReceipt,
    validation: InstanceValidation,
    encoding: serde_json::Value,
    harness_binary_blake3: String,
    solver: SolverReceipt,
    outcome: String,
    parsed: Option<ParsedSolverOutput>,
    parse_error: Option<String>,
    verification: Option<VerificationReceipt>,
    accepted_sat: bool,
    claim_scope: serde_json::Value,
}

fn hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

fn bytes_receipt(bytes: &[u8]) -> ByteReceipt {
    ByteReceipt {
        bytes: bytes.len() as u64,
        blake3: hash(bytes),
    }
}

fn read_instance(path: &Path) -> Result<(PrimeFieldSmtInstance, Vec<u8>), String> {
    let bytes = fs::read(path).map_err(|error| format!("read {}: {error}", path.display()))?;
    let instance = serde_json::from_slice(&bytes)
        .map_err(|error| format!("parse {}: {error}", path.display()))?;
    Ok((instance, bytes))
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(path)
        .map_err(|error| format!("create {} without overwrite: {error}", path.display()))?;
    file.write_all(bytes)
        .map_err(|error| format!("write {}: {error}", path.display()))
}

fn pretty<T: Serialize>(value: &T) -> Result<Vec<u8>, String> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    bytes.push(b'\n');
    Ok(bytes)
}

fn print_or_write<T: Serialize>(value: &T, out: Option<&Path>) -> Result<(), String> {
    let bytes = pretty(value)?;
    if let Some(path) = out {
        write_new(path, &bytes)
    } else {
        print!("{}", String::from_utf8_lossy(&bytes));
        Ok(())
    }
}

fn parse_for_encoding(
    output: &str,
    instance: &PrimeFieldSmtInstance,
    encoding: EncodingKind,
) -> Result<ParsedSolverOutput, String> {
    if !matches!(encoding, EncodingKind::FiniteFieldS3) {
        return parse_solver_output_for_instance(output, instance);
    }
    match parse_abscissa_output_for_instance(output, instance)? {
        ParsedAbscissaOutput::Sat { x1, x2 } => {
            let witness = lift_s3_abscissae(instance, &x1, &x2)?.ok_or_else(|| {
                format!(
                    "SAT S3 abscissae x1={x1}, x2={x2} have no independently verified signed lift"
                )
            })?;
            Ok(ParsedSolverOutput::Sat { witness })
        }
        ParsedAbscissaOutput::Unsat => Ok(ParsedSolverOutput::Unsat),
        ParsedAbscissaOutput::Unknown => Ok(ParsedSolverOutput::Unknown),
    }
}

fn solve(
    instance_path: &Path,
    solver_path: &Path,
    timeout_ms: u64,
    out_dir: &Path,
    encoding: EncodingKind,
) -> Result<bool, String> {
    let (instance, input_bytes) = read_instance(instance_path)?;
    let validation = validate_instance(&instance)?;
    let (query_text, encoding_receipt, solver_options) = match encoding {
        EncodingKind::BitVector => {
            let query = emit_smt2(&instance)?;
            (
                query.text,
                serde_json::to_value(query.receipt).map_err(|error| error.to_string())?,
                Vec::new(),
            )
        }
        EncodingKind::FiniteField => {
            let query = emit_finite_field_smt2(&instance)?;
            (
                query.text,
                serde_json::to_value(query.receipt).map_err(|error| error.to_string())?,
                vec!["--ff-solver=split".to_string()],
            )
        }
        EncodingKind::FiniteFieldS3 => {
            let query = emit_finite_field_s3_smt2(&instance)?;
            (
                query.text,
                serde_json::to_value(query.receipt).map_err(|error| error.to_string())?,
                vec!["--ff-solver=split".to_string()],
            )
        }
    };
    if let Some(parent) = out_dir.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("create {}: {error}", parent.display()))?;
    }
    fs::create_dir(out_dir).map_err(|error| {
        format!(
            "create fresh run directory {} (runs are never overwritten): {error}",
            out_dir.display()
        )
    })?;

    let frozen_input = out_dir.join("instance.json");
    let query_path = out_dir.join("query.smt2");
    let stdout_path = out_dir.join("stdout.txt");
    let stderr_path = out_dir.join("stderr.txt");
    let receipt_path = out_dir.join("receipt.json");
    write_new(&frozen_input, &input_bytes)?;
    write_new(&query_path, query_text.as_bytes())?;

    let solver_before = blake3_file(solver_path)
        .map_err(|error| format!("hash solver {}: {error}", solver_path.display()))?;
    let harness = std::env::current_exe().map_err(|error| error.to_string())?;
    let harness_binary_blake3 =
        blake3_file(&harness).map_err(|error| format!("hash harness: {error}"))?;
    let execution = run_cvc5_with_options(
        solver_path,
        &query_path,
        Duration::from_millis(timeout_ms),
        &solver_options,
    )?;
    let solver_after = blake3_file(solver_path)
        .map_err(|error| format!("rehash solver {}: {error}", solver_path.display()))?;
    write_new(&stdout_path, &execution.stdout)?;
    write_new(&stderr_path, &execution.stderr)?;

    let solver_unchanged = solver_before == solver_after;
    let mut parsed = None;
    let mut parse_error = None;
    let mut verification = None;
    let mut accepted_sat = false;
    let outcome;
    if !solver_unchanged {
        outcome = "solver_binary_changed".to_string();
    } else if execution.timed_out {
        outcome = "timeout".to_string();
    } else if !execution.success {
        outcome = "solver_error".to_string();
    } else {
        match std::str::from_utf8(&execution.stdout)
            .map_err(|error| error.to_string())
            .and_then(|output| parse_for_encoding(output, &instance, encoding))
        {
            Ok(answer) => {
                outcome = match &answer {
                    ParsedSolverOutput::Sat { witness } => {
                        let checked = verify_witness(&instance, witness)?;
                        accepted_sat = checked.verified;
                        verification = Some(checked);
                        if accepted_sat {
                            "verified_sat"
                        } else {
                            "rejected_sat_model"
                        }
                    }
                    ParsedSolverOutput::Unsat => "solver_unsat_unverified",
                    ParsedSolverOutput::Unknown => "unknown",
                }
                .to_string();
                parsed = Some(answer);
            }
            Err(error) => {
                outcome = "parse_error".to_string();
                parse_error = Some(error);
            }
        }
    }

    let receipt = SolveReceipt {
        schema: "prime-field-smt.run/v1".to_string(),
        input_path: instance_path.display().to_string(),
        input: bytes_receipt(&input_bytes),
        validation,
        encoding: encoding_receipt,
        harness_binary_blake3,
        solver: SolverReceipt {
            executable: solver_path.display().to_string(),
            blake3_before: solver_before,
            blake3_after: solver_after,
            unchanged_during_run: solver_unchanged,
            arguments: execution.arguments,
            requested_timeout_ms: timeout_ms,
            harness_grace_ms: 5_000,
            elapsed_ms: execution.elapsed.as_millis().min(u128::from(u64::MAX)) as u64,
            timed_out: execution.timed_out,
            success: execution.success,
            exit_code: execution.exit_code,
            stdout: bytes_receipt(&execution.stdout),
            stderr: bytes_receipt(&execution.stderr),
        },
        outcome: outcome.clone(),
        parsed,
        parse_error,
        verification,
        accepted_sat,
        claim_scope: json!({
            "kind": "decomposition_stage_diagnostic",
            "pollard_rho_ratio": null,
            "end_to_end_s": null,
            "speedup": null
        }),
    };
    write_new(&receipt_path, &pretty(&receipt)?)?;
    println!("{}", serde_json::to_string_pretty(&receipt).unwrap());
    Ok(matches!(
        outcome.as_str(),
        "verified_sat" | "solver_unsat_unverified"
    ))
}

fn run() -> Result<bool, String> {
    match Cli::parse().command {
        Command::Check { instance } => {
            let (instance, _) = read_instance(&instance)?;
            print_or_write(&validate_instance(&instance)?, None)?;
            Ok(true)
        }
        Command::Emit {
            instance,
            smt2,
            receipt,
            encoding,
        } => {
            let (instance, _) = read_instance(&instance)?;
            match encoding {
                EncodingKind::BitVector => {
                    let query = emit_smt2(&instance)?;
                    write_new(&smt2, query.text.as_bytes())?;
                    write_new(&receipt, &pretty(&query.receipt)?)?;
                }
                EncodingKind::FiniteField => {
                    let query = emit_finite_field_smt2(&instance)?;
                    write_new(&smt2, query.text.as_bytes())?;
                    write_new(&receipt, &pretty(&query.receipt)?)?;
                }
                EncodingKind::FiniteFieldS3 => {
                    let query = emit_finite_field_s3_smt2(&instance)?;
                    write_new(&smt2, query.text.as_bytes())?;
                    write_new(&receipt, &pretty(&query.receipt)?)?;
                }
            }
            Ok(true)
        }
        Command::Reference {
            instance,
            max_factor_points,
            out,
        } => {
            let (instance, _) = read_instance(&instance)?;
            let receipt = exhaustive_reference(&instance, max_factor_points)?;
            print_or_write(&receipt, out.as_deref())?;
            Ok(true)
        }
        Command::Verify {
            instance,
            solver_output,
            out,
            encoding,
        } => {
            let (instance, _) = read_instance(&instance)?;
            let output = fs::read_to_string(&solver_output)
                .map_err(|error| format!("read {}: {error}", solver_output.display()))?;
            let parsed = parse_for_encoding(&output, &instance, encoding)?;
            let report = match &parsed {
                ParsedSolverOutput::Sat { witness } => {
                    let verification = verify_witness(&instance, witness)?;
                    json!({"parsed": parsed, "verification": verification})
                }
                _ => json!({
                    "parsed": parsed,
                    "verification": null,
                    "note": "UNSAT/unknown has no independently checked witness"
                }),
            };
            print_or_write(&report, out.as_deref())?;
            Ok(true)
        }
        Command::Solve {
            instance,
            solver,
            timeout_ms,
            out_dir,
            encoding,
        } => solve(&instance, &solver, timeout_ms, &out_dir, encoding),
    }
}

fn main() {
    match run() {
        Ok(true) => {}
        Ok(false) => std::process::exit(2),
        Err(error) => {
            eprintln!("prime-field-smt: {error}");
            std::process::exit(1);
        }
    }
}
