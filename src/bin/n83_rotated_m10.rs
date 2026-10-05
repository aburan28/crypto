//! Command-line harness for the native n83 rotated-m10 relation gate.

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};
use crypto_lib::cryptanalysis::koblitz_rotated_chain::{
    build_n83_rotated_chain_variant, n83_public_generator, n83_public_target, write_dimacs,
    write_dimacs_model, N83CircuitVariant, N83TargetKind, OnbCurve, Type2Onb, N83, N83_RHO_SCALAR,
    N83_RHO_WALK_ITERATIONS, N83_SLOT_DIMENSIONS, N83_SUBGROUP_ORDER,
};
use serde_json::{json, Value};

const DEFAULT_NODE_CAP: usize = 8_000_000;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum TargetArg {
    Public,
    Planted,
}

impl From<TargetArg> for N83TargetKind {
    fn from(value: TargetArg) -> Self {
        match value {
            TargetArg::Public => Self::PublicRhoTarget,
            TargetArg::Planted => Self::DeterministicPlanted,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum CircuitArg {
    Complete,
    Inductive,
}

impl From<CircuitArg> for N83CircuitVariant {
    fn from(value: CircuitArg) -> Self {
        match value {
            CircuitArg::Complete => Self::CompletePerEdgeValidity,
            CircuitArg::Inductive => Self::InductiveFactorValidity,
        }
    }
}

#[derive(Debug, Parser)]
#[command(about = "Complete-point n83 rotated-m10 relation gate")]
struct Cli {
    #[command(subcommand)]
    command: Action,
}

#[derive(Debug, Subcommand)]
enum Action {
    /// Replay the exact cross-repository curve, generator, target and rho scalar.
    Check,
    /// Build the Boolean DAG and emit ordinary Tseitin DIMACS.
    Emit {
        #[arg(long, value_enum, default_value_t = TargetArg::Public)]
        target: TargetArg,
        #[arg(long, value_enum, default_value_t = CircuitArg::Complete)]
        circuit: CircuitArg,
        #[arg(long)]
        cnf: PathBuf,
        #[arg(long)]
        report: PathBuf,
        #[arg(long, default_value_t = DEFAULT_NODE_CAP)]
        node_cap: usize,
    },
    /// Emit the planted positive control plus an independently replayable model.
    CertifyPlanted {
        #[arg(long, value_enum, default_value_t = CircuitArg::Complete)]
        circuit: CircuitArg,
        #[arg(long)]
        cnf: PathBuf,
        #[arg(long)]
        model: PathBuf,
        #[arg(long)]
        report: PathBuf,
        #[arg(long, default_value_t = DEFAULT_NODE_CAP)]
        node_cap: usize,
    },
    /// Decode and independently verify a SAT model against the group law.
    Verify {
        #[arg(long, value_enum, default_value_t = TargetArg::Public)]
        target: TargetArg,
        #[arg(long, value_enum, default_value_t = CircuitArg::Complete)]
        circuit: CircuitArg,
        #[arg(long)]
        model: PathBuf,
        #[arg(long, default_value_t = DEFAULT_NODE_CAP)]
        node_cap: usize,
    },
    /// Cold-build, solve under a wall cap, and archive the complete receipt.
    Solve {
        #[arg(long, value_enum, default_value_t = TargetArg::Public)]
        target: TargetArg,
        #[arg(long, value_enum, default_value_t = CircuitArg::Complete)]
        circuit: CircuitArg,
        #[arg(long)]
        solver: PathBuf,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, default_value_t = 120)]
        seconds: u64,
        #[arg(long, default_value_t = DEFAULT_NODE_CAP)]
        node_cap: usize,
    },
}

fn write_json(path: &Path, value: &Value) -> Result<(), String> {
    let bytes = serde_json::to_vec_pretty(value).map_err(|error| error.to_string())?;
    let mut with_newline = bytes;
    with_newline.push(b'\n');
    fs::write(path, with_newline).map_err(|error| format!("write {}: {error}", path.display()))
}

fn exact_identity_report() -> Result<Value, String> {
    let field = Type2Onb::new(N83)?;
    let curve = OnbCurve::new(field);
    let generator = n83_public_generator();
    let target = n83_public_target();
    if !curve.is_on_curve(generator) || !curve.is_on_curve(target) {
        return Err("public n83 generator or target is off curve".into());
    }
    if curve.scalar_mul(generator, N83_SUBGROUP_ORDER)
        != crypto_lib::cryptanalysis::koblitz_rotated_chain::OnbPoint::INFINITY
    {
        return Err("public generator does not have the frozen subgroup order".into());
    }
    if curve.scalar_mul(generator, N83_RHO_SCALAR) != target {
        return Err("frozen rho scalar does not reproduce the public target".into());
    }
    Ok(json!({
        "schema": "n83_rotated_m10_identity_v1",
        "status": "PASS",
        "curve_id": "EC1N83Ckb1h876c2921cb64",
        "crypto_curve_id": "icv1-f2m83-tm6151469093347-debefd74",
        "model": "y^2 + x*y = x^3 + 1",
        "basis": "type_ii_optimal_normal",
        "degree": N83,
        "group_order": (4 * N83_SUBGROUP_ORDER).to_string(),
        "subgroup_order": N83_SUBGROUP_ORDER.to_string(),
        "cofactor": 4,
        "generator": generator,
        "target": target,
        "rho_recovered_scalar": N83_RHO_SCALAR.to_string(),
        "rho_walk_iterations": N83_RHO_WALK_ITERATIONS.to_string(),
        "rho_scalar_replay": true,
        "slot_dimensions": N83_SLOT_DIMENSIONS,
        "claim_boundary": "identity and representation replay only; no descent or speed claim"
    }))
}

fn chain_report(
    chain: &crypto_lib::cryptanalysis::koblitz_rotated_chain::RotatedChain,
    build_seconds: f64,
) -> Value {
    json!({
        "schema": "n83_rotated_m10_chain_v1",
        "status": "BUILT",
        "curve_id": "EC1N83Ckb1h876c2921cb64",
        "crypto_curve_id": "icv1-f2m83-tm6151469093347-debefd74",
        "target_kind": chain.target_kind,
        "circuit_variant": chain.circuit_variant,
        "target": chain.target,
        "slot_dimensions": N83_SLOT_DIMENSIONS,
        "slots": chain.slots,
        "dag": chain.counts(),
        "dag_prefix_blake3": chain.dag_prefix_blake3(),
        "build_wall_seconds": build_seconds,
        "rho_baseline": {
            "verified": true,
            "walk_iterations": N83_RHO_WALK_ITERATIONS.to_string(),
            "walk_iterations_log2": 37.55365928975891_f64,
        },
        "claim_boundary": "complete relation representation only; no natural-target yield, DLP, or rho crossover claim"
    })
}

fn build_with_report(
    target: N83TargetKind,
    circuit: N83CircuitVariant,
    node_cap: usize,
) -> Result<
    (
        crypto_lib::cryptanalysis::koblitz_rotated_chain::RotatedChain,
        Value,
    ),
    String,
> {
    let started = Instant::now();
    let chain = build_n83_rotated_chain_variant(target, circuit, node_cap)?;
    let report = chain_report(&chain, started.elapsed().as_secs_f64());
    Ok((chain, report))
}

fn main() {
    if let Err(error) = run() {
        eprintln!("n83-rotated-m10: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    match Cli::parse().command {
        Action::Check => {
            println!(
                "{}",
                serde_json::to_string_pretty(&exact_identity_report()?)
                    .map_err(|error| error.to_string())?
            );
        }
        Action::Emit {
            target,
            circuit,
            cnf,
            report,
            node_cap,
        } => {
            let (chain, mut value) = build_with_report(target.into(), circuit.into(), node_cap)?;
            let started = Instant::now();
            let receipt = write_dimacs(&chain, &cnf)?;
            value["cnf"] = serde_json::to_value(receipt).map_err(|error| error.to_string())?;
            value["cnf_write_wall_seconds"] = json!(started.elapsed().as_secs_f64());
            write_json(&report, &value)?;
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
        }
        Action::CertifyPlanted {
            circuit,
            cnf,
            model,
            report,
            node_cap,
        } => {
            let (chain, mut value) = build_with_report(
                N83TargetKind::DeterministicPlanted,
                circuit.into(),
                node_cap,
            )?;
            let cnf_receipt = write_dimacs(&chain, &cnf)?;
            let model_values = chain.known_planted_model()?;
            write_dimacs_model(&model_values, &model)?;
            let decoded = chain.decode_dimacs_model(&model)?;
            value["status"] = json!("PASS_PLANTED_CERTIFICATE");
            value["cnf"] = serde_json::to_value(cnf_receipt).map_err(|error| error.to_string())?;
            value["decoded_relation"] =
                serde_json::to_value(decoded).map_err(|error| error.to_string())?;
            write_json(&report, &value)?;
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
        }
        Action::Verify {
            target,
            circuit,
            model,
            node_cap,
        } => {
            let (chain, mut value) = build_with_report(target.into(), circuit.into(), node_cap)?;
            let decoded = chain.decode_dimacs_model(&model)?;
            value["status"] = json!("PASS_MODEL_REPLAY");
            value["decoded_relation"] =
                serde_json::to_value(decoded).map_err(|error| error.to_string())?;
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
        }
        Action::Solve {
            target,
            circuit,
            solver,
            out,
            seconds,
            node_cap,
        } => {
            if out.exists() {
                return Err(format!("refusing to overwrite {}", out.display()));
            }
            fs::create_dir_all(&out)
                .map_err(|error| format!("create {}: {error}", out.display()))?;
            let cnf = out.join("instance.cnf");
            let model = out.join("model.txt");
            let stdout_path = out.join("solver.stdout.txt");
            let stderr_path = out.join("solver.stderr.txt");
            let (chain, mut value) = build_with_report(target.into(), circuit.into(), node_cap)?;
            let cnf_receipt = write_dimacs(&chain, &cnf)?;
            value["cnf"] = serde_json::to_value(cnf_receipt).map_err(|error| error.to_string())?;
            value["solver"] = json!({
                "path": solver,
                "wall_cap_seconds": seconds,
            });
            let stdout = File::create(&stdout_path)
                .map_err(|error| format!("create {}: {error}", stdout_path.display()))?;
            let stderr = File::create(&stderr_path)
                .map_err(|error| format!("create {}: {error}", stderr_path.display()))?;
            let started = Instant::now();
            let status = Command::new(&solver)
                .arg("-t")
                .arg(seconds.to_string())
                .arg("-w")
                .arg(&model)
                .arg(&cnf)
                .stdout(Stdio::from(stdout))
                .stderr(Stdio::from(stderr))
                .status()
                .map_err(|error| format!("launch {}: {error}", solver.display()))?;
            let solver_seconds = started.elapsed().as_secs_f64();
            value["solver"]["exit_code"] = json!(status.code());
            value["solver"]["wall_seconds"] = json!(solver_seconds);
            value["status"] = match status.code() {
                Some(10) => {
                    let decoded = chain.decode_dimacs_model(&model)?;
                    value["decoded_relation"] =
                        serde_json::to_value(decoded).map_err(|error| error.to_string())?;
                    json!("SAT_VERIFIED_RELATION")
                }
                Some(20) => json!("UNSAT"),
                code => {
                    value["solver"]["nonterminal_exit"] = json!(code);
                    json!("BUDGET_OR_SOLVER_STOP")
                }
            };
            value["claim_boundary"] = json!(
                "one bounded PDP gate; no relation-rank, DLP, cold-cost crossover, or speed claim"
            );
            write_json(&out.join("result.json"), &value)?;
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
        }
    }
    Ok(())
}
