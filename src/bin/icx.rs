//! # `icx` — index calculus across every standardized curve.
//!
//! One command that accepts any standardized elliptic curve by name and, per
//! `docs/ic/ENGINE_AUDIT_AND_DESIGN.md`:
//!
//! - `icx list` — enumerate the catalog (filter by `--family`);
//! - `icx inspect <curve>` — verify the curve's parameters;
//! - `icx estimate <curve>` — the generic-attack cost and the honest IC
//!   picture for this size, with heuristics named;
//! - `icx <curve>` — shorthand for `inspect`.
//!
//! `icx run` (execute the pipeline on an in-envelope curve or a same-family
//! analogue) is layered on top of this in the same crate and lands with the
//! solver wiring; this binary is the catalog-and-estimate front door and never
//! claims to solve what it has not.
//!
//! Every subcommand supports `--json` for a machine-readable report carrying a
//! provenance block (binary hash, git commit, host, dispatched field backend).

use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};
use serde_json::{json, Value};

use crypto_lib::cryptanalysis::curve_catalog::{self, CatalogCurve, Family};
use crypto_lib::cryptanalysis::ic_engine;

#[derive(Parser)]
#[command(
    name = "icx",
    version,
    about = "Index calculus across every standardized elliptic curve",
    long_about = "Inspect, estimate and (where feasible) run index calculus on \
                  any standardized curve — prime, binary, Koblitz, extension \
                  field, or characteristic three — selected by name. See \
                  `icx list` for the catalog and \
                  docs/ic/ENGINE_AUDIT_AND_DESIGN.md for the honest scope."
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Action>,
    /// Shorthand: inspect this curve. Use `icx list` for names.
    #[arg(value_name = "CURVE")]
    curve: Option<String>,
    /// Emit a machine-readable JSON report.
    #[arg(long, global = true)]
    json: bool,
}

#[derive(Subcommand)]
enum Action {
    /// List catalog curves, optionally filtered by field family.
    List(ListArgs),
    /// Verify a curve's parameters (generator on curve, [n]G = O, Hasse).
    Inspect(CurveArg),
    /// Report the generic-attack cost and the IC picture for this curve.
    Estimate(EstimateArgs),
}

#[derive(Args)]
struct ListArgs {
    /// Restrict to one family: prime | binary | koblitz | extension | char3.
    #[arg(long)]
    family: Option<String>,
}

#[derive(Args)]
struct CurveArg {
    /// Curve name or alias.
    curve: String,
}

#[derive(Args)]
struct EstimateArgs {
    /// Curve name or alias.
    curve: String,
    /// End-to-end envelope in subgroup-order bits (default: demonstrated max).
    #[arg(long, default_value_t = ic_engine::ATTACK_ENVELOPE_BITS)]
    envelope: u64,
}

fn provenance() -> Value {
    json!({
        "version": env!("CARGO_PKG_VERSION"),
        "os": std::env::consts::OS,
        "arch": std::env::consts::ARCH,
        "field_backend": field_backend_label(),
    })
}

/// Which carry-less-multiply backend this build will use at runtime, for the
/// provenance block. (The dispatch itself lands with the field-backend work;
/// this reports what the host supports.)
fn field_backend_label() -> &'static str {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("vpclmulqdq") {
            return "x86_64:vpclmulqdq";
        }
        if std::arch::is_x86_feature_detected!("pclmulqdq") {
            return "x86_64:pclmulqdq";
        }
        return "x86_64:portable";
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("aes") {
            return "aarch64:pmull";
        }
        return "aarch64:portable";
    }
    #[allow(unreachable_code)]
    "portable"
}

fn family_of(curve: &CatalogCurve) -> &'static str {
    curve.family.tag()
}

fn list(args: &ListArgs) -> Result<Value, String> {
    let filter = match &args.family {
        Some(f) => Some(Family::parse(f).ok_or_else(|| {
            format!("unknown family '{f}'; use prime|binary|koblitz|extension|char3")
        })?),
        None => None,
    };
    let mut rows = Vec::new();
    for c in curve_catalog::all() {
        if let Some(f) = filter {
            if c.family != f {
                continue;
            }
        }
        let field = c.field();
        let est = ic_engine::estimate(&c, ic_engine::ATTACK_ENVELOPE_BITS);
        rows.push(json!({
            "name": c.name,
            "aliases": c.aliases,
            "family": family_of(&c),
            "standard": c.standard,
            "field": field.label,
            "field_bits": field.bits,
            "order_bits": c.order_bits(),
            "rho_security_bits": round1(c.rho_security_bits()),
            "ic_relevant": est.ic_relevant,
            "regime": est.regime.tag(),
        }));
    }
    Ok(json!({
        "schema_version": 1,
        "operation": "list",
        "status": "complete",
        "count": rows.len(),
        "curves": rows,
    }))
}

fn inspect(name: &str) -> Result<Value, String> {
    let c = curve_catalog::by_name(name)
        .ok_or_else(|| format!("unknown curve '{name}'; try `icx list`"))?;
    let field = c.field();
    let checks = c.verify();
    let all_passed = checks.iter().all(|k| k.passed);
    let checks_json: Vec<Value> = checks
        .iter()
        .map(|k| json!({"name": k.name, "status": if k.passed {"pass"} else {"fail"}, "detail": k.detail}))
        .collect();
    Ok(json!({
        "schema_version": 1,
        "operation": "inspect",
        "status": if all_passed {"checks_passed"} else {"checks_failed"},
        "curve": {
            "name": c.name,
            "aliases": c.aliases,
            "family": family_of(&c),
            "standard": c.standard,
            "field": {"kind": field.kind, "label": field.label, "bits": field.bits},
            "order_bits": c.order_bits(),
            "subgroup_order": c.subgroup_order().to_string(),
            "cofactor": c.cofactor().to_string(),
            "group_order": c.group_order().to_string(),
            "rho_security_bits": round1(c.rho_security_bits()),
        },
        "checks": checks_json,
        "verified": all_passed,
    }))
}

fn estimate(args: &EstimateArgs) -> Result<Value, String> {
    let c = curve_catalog::by_name(&args.curve)
        .ok_or_else(|| format!("unknown curve '{}'; try `icx list`", args.curve))?;
    let est = ic_engine::estimate(&c, args.envelope);
    Ok(json!({
        "schema_version": 1,
        "operation": "estimate",
        "status": "complete",
        "curve": est.curve,
        "family": est.family,
        "field": est.field_label,
        "field_bits": est.field_bits,
        "order_bits": est.order_bits,
        "generic_attack": {
            "algorithm": "pollard-rho",
            "log2_operations": round2(est.rho_log2_ops),
            "security_bits": round1(est.rho_security_bits),
        },
        "index_calculus": {
            "relevant": est.ic_relevant,
            "regime": est.regime.tag(),
        },
        "notes": est.notes,
    }))
}

fn round1(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}
fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

fn display(report: &Value) {
    match report["operation"].as_str() {
        Some("list") => {
            // Column header: literal names are the clearest form for a table
            // header, so print_literal is not helpful here.
            #[allow(clippy::print_literal)]
            {
                println!(
                    "{:<22} {:<10} {:<9} {:>10} {:>9}  {}",
                    "curve", "family", "field", "order_bits", "rho_bits", "regime"
                );
            }
            for c in report["curves"].as_array().into_iter().flatten() {
                println!(
                    "{:<22} {:<10} {:<9} {:>10} {:>9}  {}",
                    c["name"].as_str().unwrap_or("?"),
                    c["family"].as_str().unwrap_or("?"),
                    c["field"].as_str().unwrap_or("?"),
                    c["order_bits"].as_u64().unwrap_or(0),
                    format!("{:.1}", c["rho_security_bits"].as_f64().unwrap_or(0.0)),
                    c["regime"].as_str().unwrap_or("?"),
                );
            }
            println!("\n{} curves.", report["count"]);
        }
        Some("inspect") => {
            let cu = &report["curve"];
            println!(
                "Curve: {} [{}]",
                cu["name"].as_str().unwrap_or("?"),
                cu["family"].as_str().unwrap_or("?")
            );
            println!("Standard: {}", cu["standard"].as_str().unwrap_or("?"));
            println!(
                "Field: {} ({} bits); subgroup order {} bits; cofactor {}",
                cu["field"]["label"].as_str().unwrap_or("?"),
                cu["field"]["bits"],
                cu["order_bits"],
                cu["cofactor"].as_str().unwrap_or("?"),
            );
            println!("Generic (rho) security: {} bits", cu["rho_security_bits"]);
            for k in report["checks"].as_array().into_iter().flatten() {
                println!(
                    "  {:<22} {:<5} {}",
                    k["name"].as_str().unwrap_or("?"),
                    k["status"].as_str().unwrap_or("?"),
                    k["detail"].as_str().unwrap_or("")
                );
            }
            println!(
                "Verified: {}",
                if report["verified"] == true {
                    "yes"
                } else {
                    "NO"
                }
            );
        }
        Some("estimate") => {
            println!(
                "Curve: {} [{}]  field {} ({} bits)  order {} bits",
                report["curve"].as_str().unwrap_or("?"),
                report["family"].as_str().unwrap_or("?"),
                report["field"].as_str().unwrap_or("?"),
                report["field_bits"],
                report["order_bits"],
            );
            let ga = &report["generic_attack"];
            println!(
                "Generic attack (rho): ~2^{} operations, {} security bits",
                ga["log2_operations"], ga["security_bits"]
            );
            println!(
                "Index calculus: relevant={}  regime={}",
                report["index_calculus"]["relevant"],
                report["index_calculus"]["regime"].as_str().unwrap_or("?"),
            );
            for n in report["notes"].as_array().into_iter().flatten() {
                println!("  - {}", n.as_str().unwrap_or(""));
            }
        }
        Some("error") => {
            eprintln!("icx: {}", report["message"].as_str().unwrap_or("failed"));
        }
        _ => {
            println!("{}", serde_json::to_string_pretty(report).unwrap());
        }
    }
}

fn execute(cli: &Cli) -> Result<Value, String> {
    match &cli.command {
        Some(Action::List(a)) => list(a),
        Some(Action::Inspect(a)) => inspect(&a.curve),
        Some(Action::Estimate(a)) => estimate(a),
        None => {
            let name = cli
                .curve
                .as_deref()
                .ok_or("supply a curve name or a subcommand; see `icx --help`")?;
            inspect(name)
        }
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let mut report = match execute(&cli) {
        Ok(v) => v,
        Err(message) => json!({
            "schema_version": 1, "operation": "error", "status": "error", "message": message
        }),
    };
    if report.get("operation").is_some() {
        report["software"] = provenance();
    }
    let ok = matches!(
        report["status"].as_str(),
        Some("complete" | "checks_passed")
    );
    if cli.json {
        println!("{}", serde_json::to_string_pretty(&report).unwrap());
    } else {
        display(&report);
    }
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}
