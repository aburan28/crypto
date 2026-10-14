//! `ic prime`: the exact pair-sum yield of a small-x factor base on a
//! prime-field curve.
//!
//! Nothing else in `ic` selects a factor base over a prime field: `search`,
//! `run` and `workflow` are binary Koblitz pipelines.  This command measures
//! the one number a prime-field index calculus with two-point sums turns on,
//! the probability that a random group element is a signed sum of two base
//! points, exactly (from every pair sum) and by Monte Carlo where the yield is
//! large enough to see.  It runs on NIST P-192 and on toy primes, where the
//! same count can be checked against brute force and against hits.
use clap::{Args, ValueEnum};
use crypto_lib::cryptanalysis::{
    binary_small_x::{self, BinCurve},
    prime_small_x::{measure, toy_curve, Curve, Options},
    secg_standards::{secg, SecgField},
};
use serde::Serialize;
use serde_json::{json, Value};
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PrimeBase {
    /// The points with the smallest abscissae `x = 0, 1, 2, …` on the curve.
    SmallX,
}

#[derive(Args, Clone, Debug, Serialize)]
pub struct PrimeArgs {
    /// A SECG curve of at most 192 bits over a prime field or of odd degree
    /// at most 256 over a binary field (`secp192r1` also as `p192`,
    /// `nist-p192`; `sect113r1`, `sect163r2`, `secp160k1`, …), or a toy
    /// curve: `toy:<bits>[,<seed>]` over a prime just below `2^bits` (8 to
    /// 28 bits), `toy-binary:<degree>[,<seed>]` over `F_2^degree` (11, 13,
    /// 17 or 19).
    #[arg(long, default_value = "secp192r1")]
    pub curve: String,
    /// Factor-base family; the small-x base is the only one so far.
    #[arg(long, value_enum, default_value_t = PrimeBase::SmallX)]
    pub factor_base: PrimeBase,
    /// Abscissae in the base (the base has twice as many points). The pair
    /// sums are about `size²` and are held in memory.
    #[arg(long, default_value_t = 2048, value_parser = clap::value_parser!(u32).range(2..=8192))]
    pub size: u32,
    /// Random group elements tested by Monte Carlo.
    #[arg(long, default_value_t = 1_000_000, value_parser = clap::value_parser!(u64).range(0..=1_000_000_000))]
    pub trials: u64,
    /// Planted decompositions that must be detected and recovered.
    #[arg(long, default_value_t = 32, value_parser = clap::value_parser!(u32).range(0..=4096))]
    pub controls: u32,
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
}

pub fn run(args: PrimeArgs) -> Result<Value, String> {
    let started = Instant::now();
    let opts = Options {
        size: args.size as usize,
        trials: args.trials,
        controls: args.controls as usize,
        seed: args.seed,
    };
    let name = args.curve.to_ascii_lowercase();
    let toy_args = |spec: &str, what: &str| -> Result<(u32, u64), String> {
        let mut parts = spec.split(',');
        let n: u32 = parts
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| format!("{what} curves are written {what}:<size>[,<seed>]"))?;
        let seed: u64 = parts.next().map_or(Ok(args.seed), |s| {
            s.parse().map_err(|_| "toy curve seed must be an integer")
        })?;
        Ok((n, seed))
    };
    let report = if matches!(name.as_str(), "secp192r1" | "p192" | "nist-p192" | "p-192") {
        measure(&Curve::p192(), &opts)?
    } else if let Some(spec) = name.strip_prefix("toy-binary:") {
        let (degree, seed) = toy_args(spec, "toy-binary")?;
        binary_small_x::measure(&binary_small_x::toy_curve(degree, seed)?, &opts)?
    } else if let Some(spec) = name.strip_prefix("toy:") {
        let (bits, seed) = toy_args(spec, "toy")?;
        measure(&toy_curve(bits, seed)?, &opts)?
    } else if let Some(c) = secg(&name) {
        match &c.field {
            SecgField::Binary { .. } => binary_small_x::measure(&BinCurve::from_secg(c)?, &opts)?,
            SecgField::Prime { .. } => measure(&Curve::from_secg(c)?, &opts)?,
        }
    } else {
        return Err(format!(
            "unknown curve {:?}; use a SECG name from `ic list` (primes up to 192 bits, binary fields of odd degree up to 256), secp192r1, toy:<bits> or toy-binary:<degree>",
            args.curve
        ));
    };
    Ok(json!({
        "schema_version": 1,
        "operation": "prime",
        "status": "complete",
        "evidence_scope": "exact_pair_sum_yield_with_monte_carlo",
        "arguments": args,
        "report": report,
        "elapsed_seconds": started.elapsed().as_secs_f64(),
        "limitations": [
            "Two-point sums only: a relation is a signed sum of two base points of the small-x base.",
            "The yield is exact for the base built. The Monte Carlo rate agrees with it only where the yield is large enough to produce hits; at P-192 it is a null result and bounds the rate from above.",
            "The cost of each trial, deciding by a Semaev S₃ solve whether a target decomposes, is not part of this measurement, and neither is any linear algebra or descent: a relation yield is not a discrete logarithm.",
            "The base is capped at 8192 abscissae by the memory of the pair-sum set.",
            "With a cofactor h > 1 (every SECG binary curve, secp112r2, secp128r2) each small-x point enters the base as [h]P and targets are [h] times a uniform point, so base and targets lie in the prime-order subgroup and the yield is over that subgroup; keeping raw small-x points and counting only the pairs that sum into it would give about 1/h of this yield."
        ]
    }))
}
