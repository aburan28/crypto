//! `ic boundary`: the index-calculus boundary ledger.
//!
//! Runs the three regimes of
//! [`crypto_lib::cryptanalysis::ic_boundary`] — generic prime field,
//! generic binary field, Koblitz — over their size ladders, prices every
//! phase of every variant in one unit against the generic floor and a
//! counted Pollard rho, fits the exponents, and writes it all as one
//! JSON report (with the Markdown table inside it).
//!
//!     ic boundary --quick --json
//!     ic boundary --regime koblitz --koblitz-degrees 23,31 --repeats 2 --out ledger.json
//!     ic boundary --out docs/ic/runs/ic-boundary-ledger-YYYY-MM-DD.json

use clap::{Args, ValueEnum};
use crypto_lib::cryptanalysis::ic_boundary::{
    fit_exponents, format_markdown, run_char2_ladder, run_koblitz_ladder, run_prime_ladder,
    BoundaryConfig, BoundaryLedger,
};
use crypto_lib::cryptanalysis::ic_oracle_pricing::{
    format_oracle_markdown, price_oracles, OraclePricingConfig,
};
use serde_json::{json, Value};
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Regime {
    All,
    Prime,
    Char2,
    Koblitz,
}

#[derive(Args, Clone, Debug)]
pub struct BoundaryArgs {
    /// Which regime to run.
    #[arg(long, value_enum, default_value_t = Regime::All)]
    pub regime: Regime,
    /// Prime-field ladder, in subgroup bits (roster curves where the
    /// roster has one, generated prime-order curves otherwise).
    #[arg(long, value_delimiter = ',')]
    pub prime_bits: Option<Vec<u32>>,
    /// Generic binary ladder, in field degrees (random curves).
    #[arg(long, value_delimiter = ',')]
    pub char2_degrees: Option<Vec<u32>>,
    /// Koblitz ladder, in field degrees.
    #[arg(long, value_delimiter = ',')]
    pub koblitz_degrees: Option<Vec<u32>>,
    /// Repeats per instance; rho and every variant see the same target
    /// within a repeat.
    #[arg(long)]
    pub repeats: Option<usize>,
    #[arg(long)]
    pub seed: Option<u64>,
    /// The small configuration the tests use.
    #[arg(long)]
    pub quick: bool,
    /// Largest Koblitz degree at which the un-folded control runs.
    #[arg(long)]
    pub no_fold_max_degree: Option<u32>,
    /// Largest degree at which the S₄ pairs-and-solve oracle runs.
    #[arg(long)]
    pub s4_max_degree: Option<u32>,
    /// Also price the decomposition oracles (enumerate, meet in the
    /// middle, S₄, matrix-F4, SAT) per target on the Koblitz Semaev
    /// systems, with the first fall degree of each system.
    #[arg(long)]
    pub oracles: bool,
    /// Oracle cells as `n:m` pairs, e.g. `9:2,15:3`; default ladder when absent.
    #[arg(long, value_delimiter = ',')]
    pub oracle_cells: Option<Vec<String>>,
    /// Targets per oracle cell.
    #[arg(long)]
    pub oracle_targets: Option<usize>,
    /// Run the walk rows with one jump table for every segment and the
    /// guard cleared at each restart, the walk as the Round-2 ladder first
    /// ran it: segments merge like rho's walks and a repeated
    /// decomposable target pins the logarithm.  For the merge diagnostic
    /// only; the default draws fresh jumps per segment.
    #[arg(long)]
    pub walk_shared_jumps: bool,
}

fn host() -> Value {
    let cpu = std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("model name"))
                .and_then(|l| l.split(':').nth(1))
                .map(|v| v.trim().to_string())
        });
    let commit = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
    json!({
        "cpu": cpu,
        "logical_cores": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0),
        "os": std::env::consts::OS,
        "arch": std::env::consts::ARCH,
        "git_commit": commit,
    })
}

pub fn run(args: BoundaryArgs, json: bool) -> Result<Value, String> {
    let mut cfg = if args.quick {
        BoundaryConfig::quick()
    } else {
        BoundaryConfig::default()
    };
    if let Some(v) = args.prime_bits {
        cfg.prime_bits = v;
    }
    if let Some(v) = args.char2_degrees {
        cfg.char2_degrees = v;
    }
    if let Some(v) = args.koblitz_degrees {
        cfg.koblitz_degrees = v;
    }
    if let Some(v) = args.repeats {
        cfg.repeats = v.max(1);
    }
    if let Some(v) = args.seed {
        cfg.seed = v;
    }
    if let Some(v) = args.no_fold_max_degree {
        cfg.koblitz_no_fold_max_degree = v;
    }
    if let Some(v) = args.s4_max_degree {
        cfg.s4_max_degree = v;
    }
    cfg.walk_shared_jumps = args.walk_shared_jumps;
    for &b in &cfg.prime_bits {
        if !(8..=32).contains(&b) {
            return Err(format!("prime ladder bits must lie in 8..=32, got {b}"));
        }
    }
    for &n in cfg.char2_degrees.iter().chain(&cfg.koblitz_degrees) {
        if !(5..=62).contains(&n) {
            return Err(format!("field degrees must lie in 5..=62, got {n}"));
        }
    }

    let started = Instant::now();
    let mut progress = |line: &str| {
        if !json {
            eprintln!("  {line}");
        }
    };
    let mut instances = Vec::new();
    if matches!(args.regime, Regime::All | Regime::Prime) {
        instances.extend(run_prime_ladder(&cfg, &mut progress));
    }
    if matches!(args.regime, Regime::All | Regime::Char2) {
        instances.extend(run_char2_ladder(&cfg, &mut progress));
    }
    if matches!(args.regime, Regime::All | Regime::Koblitz) {
        instances.extend(run_koblitz_ladder(&cfg, &mut progress));
    }
    let oracle_pricing = if args.oracles {
        let mut ocfg = if args.quick {
            OraclePricingConfig::quick()
        } else {
            OraclePricingConfig::default()
        };
        ocfg.seed = cfg.seed;
        if let Some(cells) = &args.oracle_cells {
            let mut parsed = Vec::new();
            for c in cells {
                let (n, m) = c
                    .split_once(':')
                    .ok_or_else(|| format!("oracle cell `{c}` is not of the form n:m"))?;
                let n: u32 = n.parse().map_err(|_| format!("bad degree in `{c}`"))?;
                let m: u32 = m.parse().map_err(|_| format!("bad summand count in `{c}`"))?;
                if !(5..=62).contains(&n) || !(2..=4).contains(&m) {
                    return Err(format!("oracle cell `{c}` out of range"));
                }
                parsed.push((n, m));
            }
            ocfg.cells = parsed;
        }
        if let Some(t) = args.oracle_targets {
            ocfg.targets = t.max(1);
        }
        let cells = price_oracles(&ocfg, &mut progress);
        let all_agree = cells.iter().all(|c| c.disagreements == 0);
        Some(json!({
            "config": ocfg,
            "all_agree": all_agree,
            "cells": cells,
            "markdown": format_oracle_markdown(&cells),
        }))
    } else {
        None
    };
    let fits = fit_exponents(&instances);
    let ledger = BoundaryLedger {
        unit: "S = group-addition equivalents / sqrt(r); native counts exact, conversions measured per instance".into(),
        floor: "generic: sqrt(pi/(2A)) in S, A = usable automorphisms; relations: trials >= (K+1)/min(1, C(F+m-1,m)/#E)".into(),
        reference: "Pollard rho on the same instance, same process, exact operation count, verified [d]G = Q".into(),
        instances,
        fits,
    };
    let all_verified = ledger
        .instances
        .iter()
        .all(|i| i.rho_verified_all && i.variants.iter().all(|v| v.verified));
    let oracles_agree = oracle_pricing
        .as_ref()
        .map_or(true, |o| o["all_agree"] == true);
    let markdown = format_markdown(&ledger);
    let ran_something = !ledger.instances.is_empty()
        || oracle_pricing
            .as_ref()
            .is_some_and(|o| o["cells"].as_array().is_some_and(|c| !c.is_empty()));
    let status = if ran_something && all_verified && oracles_agree {
        "complete"
    } else {
        "incomplete"
    };
    Ok(json!({
        "schema_version": 1,
        "operation": "boundary",
        "status": status,
        "what_this_is": "Every index-calculus variant of three regimes priced end to end in one unit (group-addition equivalents per sqrt of the subgroup order) against the generic floor and a counted Pollard rho on the same instance; native counts are exact, conversion factors are measured on the host and recorded.",
        "what_this_is_not": [
            "not a claim about any deployed curve: the largest instance is a 40-bit Koblitz subgroup",
            "not a wall-clock benchmark: wall times are carried as a practicality note, the metric is the operation count",
            "not a best-of-breed composite: every row is one process on one instance with every phase inside it",
            "extrapolations, where a reader makes them from the fitted exponents, are extrapolations"
        ],
        "config": cfg,
        "host": host(),
        "elapsed_seconds": started.elapsed().as_secs_f64(),
        "all_verified": all_verified,
        "ledger": ledger,
        "markdown": markdown,
        "oracle_pricing": oracle_pricing,
    }))
}
