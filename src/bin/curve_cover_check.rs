//! Construct and verify same-field hyperelliptic covers of catalog models.
#[path = "curve_cover_check/checker.rs"]
mod checker;
#[path = "curve_cover_check/links.rs"]
mod links;
#[path = "curve_cover_check/models.rs"]
mod models;
#[cfg(test)]
#[path = "curve_cover_check/tests.rs"]
mod tests;

use clap::Parser;
use std::{fs, path::PathBuf};

#[derive(Parser)]
#[command(about = "Construct and replay explicit H -> E cover certificates")]
struct Args {
    #[arg(long, default_value = "docs/curves/registry.json")]
    registry: PathBuf,
    #[arg(long, default_value = "docs/curves/covers.json")]
    output: PathBuf,
    /// Recompute the entire report and require byte-identical saved output.
    #[arg(long)]
    check: bool,
    /// Write the HC1/CV1/ICV1/EC1 graph as YAML 1.2 (JSON-compatible syntax).
    #[arg(long, default_value = "docs/curves/cover-links.yaml")]
    links: PathBuf,
}

fn run() -> Result<(), String> {
    let args = Args::parse();
    if args.registry == args.output
        || matches!((args.registry.canonicalize(), args.output.canonicalize()), (Ok(a), Ok(b)) if a == b)
    {
        return Err("registry and output must be different files".into());
    }
    let input = fs::read(&args.registry).map_err(|e| e.to_string())?;
    let report = checker::catalog(&input)?;
    let registry: serde_json::Value = serde_json::from_slice(&input).map_err(|e| e.to_string())?;
    let graph = links::graph(&registry, &report)?;
    let graph_text = serde_json::to_string_pretty(&graph).map_err(|e| e.to_string())? + "\n";
    if [&args.registry, &args.output].iter().any(|p| {
        *p == &args.links
            || matches!((p.canonicalize(), args.links.canonicalize()), (Ok(a), Ok(b)) if a == b)
    }) {
        return Err("registry, report and graph paths must be distinct".into());
    }
    if args.check {
        if fs::read_to_string(&args.links).map_err(|e| e.to_string())? != graph_text {
            return Err("cover linkage graph is stale or altered".into());
        }
    } else {
        fs::write(&args.links, graph_text).map_err(|e| e.to_string())?;
    }
    let text = serde_json::to_string_pretty(&report).map_err(|e| e.to_string())? + "\n";
    if args.check {
        let saved = fs::read_to_string(&args.output).map_err(|e| e.to_string())?;
        if saved != text {
            return Err("cover report is stale or altered; regenerate and review it".into());
        }
    } else {
        fs::write(&args.output, text).map_err(|e| e.to_string())?;
    }
    println!("{}", report["summary"]);
    if report["summary"]["invalid_input"].as_u64().unwrap_or(0) > 0 {
        return Err("invalid catalog models; findings retained for inspection".into());
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("curve_cover_check: {e}");
        std::process::exit(1);
    }
}
