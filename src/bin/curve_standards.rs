//! Offline, source-pinned standard-parameter importer. All curve arithmetic,
//! model identities and generator checks run natively.
#[allow(dead_code)]
#[path = "curve_cover_check/checker.rs"]
mod checker;
#[path = "curve_standards/import.rs"]
mod import;
#[path = "curve_cover_check/models.rs"]
mod models;
use clap::Parser;
use serde_json::{json, Value};
use std::{fs, path::PathBuf};
#[derive(Parser)]
struct Args {
    #[arg(long, default_value = "docs/curves/standards")]
    directory: PathBuf,
    #[arg(long)]
    check: bool,
}
fn run() -> Result<(), String> {
    let args = Args::parse();
    let mut inputs = vec![];
    let mut source_hashes = serde_json::Map::new();
    for name in ["parameters.json", "supplemental.json"] {
        let bytes = fs::read(args.directory.join(name)).map_err(|e| e.to_string())?;
        source_hashes.insert(name.into(), json!(checker::digest(&bytes)));
        let doc: Value = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        inputs.extend(
            doc["curves"]
                .as_array()
                .ok_or("missing parameters")?
                .iter()
                .cloned(),
        );
    }
    let (rows, findings) = import::catalog(&inputs)?;
    let imported = findings
        .iter()
        .filter(|r| r["status"] == "imported")
        .count();
    let coverage = json!({"schema_version":"standards-coverage/v1","source_sha256":source_hashes,
        "global_completeness":"not_established","scope":"Every row in the pinned source inventories; inclusion is not an assertion of agency adoption or security approval",
        "summary":{"source_records":inputs.len(),"imported":imported,"unresolved":inputs.len()-imported},"records":findings});
    for (name, doc) in [
        (
            "registry.json",
            json!({"schema_version":"standard-registry/v1","source_sha256":source_hashes,"curves":rows}),
        ),
        ("coverage.json", coverage.clone()),
    ] {
        let text = serde_json::to_string_pretty(&doc).map_err(|e| e.to_string())? + "\n";
        let path = args.directory.join(name);
        if args.check {
            if fs::read_to_string(path).map_err(|e| e.to_string())? != text {
                return Err(format!("stale standards output: {name}"));
            }
        } else {
            fs::write(path, text).map_err(|e| e.to_string())?;
        }
    }
    println!("{}", coverage["summary"]);
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("curve_standards: {e}");
        std::process::exit(1);
    }
}
