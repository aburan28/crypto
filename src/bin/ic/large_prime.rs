use std::path::PathBuf;

use clap::Args;
use crypto_lib::cryptanalysis::ecbench_large_prime::{
    import_family, load_manifest_with_source, solve, SolveConfig,
};
use serde_json::{json, Value};

#[derive(Clone, Debug, Args)]
pub struct LargePrimeArgs {
    /// ECBench curve/point corpus manifest.
    #[arg(long)]
    pub manifest: PathBuf,
    /// Family id to solve; defaults to the first family after validating all.
    #[arg(long)]
    pub family: Option<String>,
    /// Validate the full curve/point battery without running the solver.
    #[arg(long)]
    pub validate_only: bool,
    /// Exact relation arity: an integer or n-1.
    #[arg(long, default_value = "n-1")]
    pub summands: String,
    /// Maximum accepted large primes per partial relation (0, 1, or 2).
    #[arg(long, default_value_t = 2)]
    pub large_primes: u8,
    /// Override the corpus small factor-base dimension.
    #[arg(long)]
    pub small_dimension: Option<u32>,
    /// Override the corpus factor-base envelope dimension.
    #[arg(long)]
    pub envelope_dimension: Option<u32>,
    /// Relation-trial cap.
    #[arg(long, default_value_t = 100_000)]
    pub max_trials: u64,
    /// Exact meet-in-the-middle state cap.
    #[arg(long, default_value_t = 1_000_000)]
    pub max_states: u64,
    /// Deterministic relation-search seed.
    #[arg(long, default_value_t = 7)]
    pub seed: u64,
}

pub fn run(args: LargePrimeArgs) -> Result<Value, String> {
    let loaded = load_manifest_with_source(&args.manifest)?;
    let manifest_sha256 = loaded.sha256;
    let manifest_bytes = loaded.bytes;
    let manifest_schema = loaded.manifest.schema.clone();
    let manifest = loaded.manifest;
    let wanted = args
        .family
        .clone()
        .unwrap_or_else(|| manifest.families[0].id.clone());
    let mut battery = Vec::with_capacity(manifest.families.len());
    let mut selected = None;
    for family in &manifest.families {
        let (imported, facts) = import_family(family)?;
        battery.push(serde_json::to_value(facts).map_err(|e| e.to_string())?);
        if family.id == wanted {
            selected = Some(imported);
        }
    }
    let imported = selected.ok_or_else(|| {
        let ids = manifest
            .families
            .iter()
            .map(|f| f.id.as_str())
            .collect::<Vec<_>>()
            .join(", ");
        format!("unknown family `{wanted}`; available: {ids}")
    })?;
    if args.validate_only {
        return Ok(json!({
            "schema_version": 1,
            "operation": "large-prime",
            "status": "checks_passed",
            "manifest": args.manifest,
            "manifest_sha256": manifest_sha256,
            "manifest_bytes": manifest_bytes,
            "manifest_schema": manifest_schema,
            "battery": battery,
            "run": Value::Null,
        }));
    }
    let summands = match args.summands.as_str() {
        "n-1" => imported.instance.n - 1,
        value => value
            .parse::<u32>()
            .map_err(|_| format!("--summands must be an integer or n-1, not `{value}`"))?,
    };
    let config = SolveConfig {
        small_dimension: args.small_dimension.unwrap_or(imported.small_dimension),
        envelope_dimension: args
            .envelope_dimension
            .unwrap_or(imported.envelope_dimension),
        summands,
        max_large_primes: args.large_primes,
        max_trials: args.max_trials,
        max_states: args.max_states,
        seed: args.seed,
    };
    // `expected` is deliberately read only after the target-only solver call.
    let report = solve(&imported.instance, imported.target, config)?;
    let expected = imported.expected;
    if report.recovered.is_some_and(|d| d != expected) {
        return Err(format!(
            "solver returned {}, but the post-solve corpus check expected {expected}",
            report.recovered.unwrap()
        ));
    }
    let status = if report.recovered == Some(expected)
        && report.verified_fast
        && report.verified_independent
    {
        "complete"
    } else {
        "stopped"
    };
    Ok(json!({
        "schema_version": 1,
        "operation": "large-prime",
        "status": status,
        "manifest": args.manifest,
        "manifest_sha256": manifest_sha256,
        "manifest_bytes": manifest_bytes,
        "manifest_schema": manifest_schema,
        "battery": battery,
        "run": {
            "family": imported.id,
            "expected": expected,
            "matches_expected": report.recovered == Some(expected),
            "report": report,
        },
        "scope": "correctness experiment; native combinatorics and modular elimination are counted but unpriced, so this is not a performance or scaling claim",
    }))
}
