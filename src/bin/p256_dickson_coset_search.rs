//! Frozen eight-coset P-256 Dickson factor-base search.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, RelationMetrics, WideBuildResult, CURVE_SLUG,
};
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::One;
use serde::Serialize;
use serde_json::json;

const DOMAIN: &str =
    "icv1-fp256-t89188191154553853111372247798585809583-f188c491/dickson-terminal/v1/";
const DEPTH: u32 = 18;
const CANDIDATES: u32 = 8;

#[derive(Parser)]
#[command(about = "Search the frozen hash-derived P-256 Dickson cosets")]
struct Cli {
    /// Registered ICV1 slug.  The exact P-256 model is the only accepted curve.
    #[arg(long)]
    curve: String,
    /// Compact result JSON path; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Rebuild the selected fibre and compare every emitted point row.
    #[arg(long)]
    verify: bool,
}

#[derive(Clone, Debug, Serialize)]
struct Candidate {
    index: u32,
    root_exponent: String,
    terminal: String,
    fb_id: String,
    fb_sha256: String,
    columns: u64,
    signed_points: u64,
    points_sha256: String,
}

#[derive(Clone, Debug, Serialize)]
struct SearchResult {
    schema: String,
    curve: String,
    family: String,
    depth: u32,
    derivation_domain: String,
    candidates: Vec<Candidate>,
    selected_index: u32,
    selected_fb1_preimage: serde_json::Value,
    selected_relation_m17: RelationMetrics,
    verified: bool,
}

fn root_exponent(index: u32) -> BigUint {
    let digest = sha256(format!("{DOMAIN}{index}").as_bytes());
    let mask = (BigUint::one() << 78usize) - BigUint::one();
    (BigUint::from_bytes_be(&digest) & mask) | BigUint::one()
}

fn candidate(index: u32, result: &WideBuildResult) -> Result<Candidate, String> {
    let fb = &result.dump.factor_base;
    Ok(Candidate {
        index,
        root_exponent: fb
            .params
            .get("root_exponent")
            .ok_or("Dickson result has no root exponent".to_string())?
            .clone(),
        terminal: fb
            .params
            .get("terminal")
            .ok_or("Dickson result has no terminal".to_string())?
            .clone(),
        fb_id: fb.fb_id.clone(),
        fb_sha256: fb.fb_sha256.clone(),
        columns: fb.columns,
        signed_points: fb.signed_points,
        points_sha256: fb.points_sha256.clone(),
    })
}

fn run(cli: Cli) -> Result<(), String> {
    if cli.curve != CURVE_SLUG {
        return Err(format!(
            "P-256 Dickson coset search has no registered curve `{}`",
            cli.curve
        ));
    }
    let mut candidates = Vec::with_capacity(CANDIDATES as usize);
    let mut selected: Option<(u32, WideBuildResult)> = None;
    for index in 0..CANDIDATES {
        let exponent = root_exponent(index);
        let spec = format!(
            "dickson-torus:depth={DEPTH},root_exponent=0x{}",
            exponent.to_str_radix(16)
        );
        let result = p256_dickson_factor_base::build(&spec)?;
        let summary = candidate(index, &result)?;
        if selected.as_ref().is_none_or(|(best_index, best)| {
            summary.columns > best.dump.factor_base.columns
                || (summary.columns == best.dump.factor_base.columns && index < *best_index)
        }) {
            selected = Some((index, result));
        }
        candidates.push(summary);
    }
    let (selected_index, selected) = selected.ok_or("no P-256 Dickson candidates".to_string())?;
    if cli.verify {
        p256_dickson_factor_base::verify(&selected.dump)?;
    }
    let fb = &selected.dump.factor_base;
    let preimage = json!({
        "schema": p256_dickson_factor_base::FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": p256_dickson_factor_base::FAMILY,
        "params": fb.params,
        "columns": fb.columns,
        "signed_points": fb.signed_points,
        "points_sha256": fb.points_sha256,
    });
    let result = SearchResult {
        schema: "p256.dickson_coset_search/v1".into(),
        curve: CURVE_SLUG.into(),
        family: p256_dickson_factor_base::FAMILY.into(),
        depth: DEPTH,
        derivation_domain: DOMAIN.into(),
        candidates,
        selected_index,
        selected_fb1_preimage: preimage,
        selected_relation_m17: p256_dickson_factor_base::relation_metrics(fb.columns, 17)?,
        verified: cli.verify,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match cli.out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    eprintln!(
        "{}: selected k={} with {} columns{}",
        CURVE_SLUG,
        result.selected_index,
        result.selected_relation_m17.columns,
        if result.verified { " (verified)" } else { "" },
    );
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
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
    fn exponents_are_odd_bounded_and_distinct() {
        let limit = BigUint::one() << 78usize;
        let values: Vec<BigUint> = (0..CANDIDATES).map(root_exponent).collect();
        assert!(values.iter().all(|value| value < &limit && value.bit(0)));
        for pair in values.iter().enumerate() {
            assert!(!values[..pair.0].contains(pair.1));
        }
    }
}
