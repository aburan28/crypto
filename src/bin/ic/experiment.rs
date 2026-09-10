//! Bounded experiments on internally generated, known-answer toy instances.
use super::params::{self, Field, Fixture, Parameters};
use clap::{Args, ValueEnum};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    factor_x_n_minus_1, koblitz_index_calculus_dlp_with_progress, order_of_2_mod_n,
    DecompositionStrategy, KoblitzCurve, KoblitzIcEvent, KoblitzIcOptions, MAX_N,
};
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Serialize;
use serde_json::{json, Value};
use std::{
    io::{Read, Write},
    process::{Command, Stdio},
    time::{Duration, Instant},
};

pub const MAX_FACTOR_DIMENSION: u32 = 12;
pub fn degree(value: &str) -> Result<u32, String> {
    let n = value
        .parse::<u32>()
        .map_err(|_| "degree must be an integer")?;
    if n < 3 || n > MAX_N || n % 2 == 0 {
        return Err(format!(
            "synthetic degree must be odd and 3..={}",
            MAX_N - 1
        ));
    }
    Ok(n)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Solver {
    Groebner,
    Sat,
    Enumerate,
}
impl Solver {
    pub fn name(self) -> &'static str {
        match self {
            Self::Groebner => "groebner",
            Self::Sat => "sat",
            Self::Enumerate => "enumerate",
        }
    }
    fn strategy(self) -> DecompositionStrategy {
        match self {
            Self::Groebner => DecompositionStrategy::Groebner,
            Self::Sat => DecompositionStrategy::Sat,
            Self::Enumerate => DecompositionStrategy::Enumerate,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Args, Serialize)]
pub struct RunArgs {
    /// Odd extension degree for a generated GF(2^n) test curve, 3..=23.
    #[arg(long,default_value_t=9,value_parser=degree)]
    pub degree: u32,
    /// Koblitz coefficient a; b is always 1.
    #[arg(long,default_value_t=0,value_parser=clap::value_parser!(u8).range(0..=1))]
    pub curve_a: u8,
    /// Public logarithm; defaults to 53 unless --random-target is selected.
    #[arg(long,value_parser=clap::value_parser!(u64).range(1..),conflicts_with="random_target")]
    pub known_log: Option<u64>,
    /// Generate a deterministic random known-answer target from --seed.
    #[arg(long)]
    pub random_target: bool,
    #[arg(long, default_value_t = 0x4b_6f_62_6c_69_74_7a_00u64)]
    pub seed: u64,
    /// Candidate index in the existing degree-ord_n(2) factor-base family.
    #[arg(long, default_value_t = 0)]
    pub factor_index: usize,
    #[arg(long,default_value_t=20_000,value_parser=clap::value_parser!(u32).range(1..=20_000))]
    pub max_trials: u32,
    #[arg(long,value_enum,default_value_t=Solver::Groebner)]
    pub solver: Solver,
}
impl Default for RunArgs {
    fn default() -> Self {
        Self {
            degree: 9,
            curve_a: 0,
            known_log: None,
            random_target: false,
            seed: 0x4b_6f_62_6c_69_74_7a_00,
            factor_index: 0,
            max_trials: 20_000,
            solver: Solver::Groebner,
        }
    }
}
#[derive(Clone, Debug, Args)]
pub struct GenerateArgs {
    #[arg(long,default_value_t=9,value_parser=degree)]
    pub degree: u32,
    #[arg(long,default_value_t=0,value_parser=clap::value_parser!(u8).range(0..=1))]
    pub curve_a: u8,
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
}
#[derive(Clone, Debug, Args)]
pub struct CompareArgs {
    #[arg(long,default_value_t=7,value_parser=degree)]
    pub degree: u32,
    #[arg(long,default_value_t=1,value_parser=clap::value_parser!(u8).range(0..=1))]
    pub curve_a: u8,
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
    /// Identical generated training fixtures per candidate.
    #[arg(long,default_value_t=3,value_parser=clap::value_parser!(u32).range(1..=8))]
    pub samples: u32,
    /// Independent generated fixtures used to check the training winner.
    #[arg(long,default_value_t=2,value_parser=clap::value_parser!(u32).range(1..=8))]
    pub holdout: u32,
    #[arg(long,default_value_t=2000,value_parser=clap::value_parser!(u32).range(1..=20_000))]
    pub max_trials: u32,
    #[arg(long,default_value_t=30,value_parser=clap::value_parser!(u32).range(1..=120))]
    pub timeout_seconds: u32,
    #[arg(long,value_enum,default_value_t=Solver::Enumerate)]
    pub solver: Solver,
}
fn curve(n: u32, a: u8) -> Result<KoblitzCurve, String> {
    // Keep this guard even for internal callers, independently of Clap.
    if n < 3 || n > MAX_N || n % 2 == 0 || a > 1 {
        return Err("unsupported synthetic curve parameters".into());
    }
    KoblitzCurve::new(a, n).ok_or_else(|| {
        format!("K_{a} / GF(2^{n}) has no usable prime-order subgroup in the existing constructor")
    })
}
fn known(curve: &KoblitzCurve, args: &RunArgs) -> Result<BigUint, String> {
    let n = if args.random_target {
        let mut rng = StdRng::seed_from_u64(args.seed ^ 0x534f_4c56_4552_5447);
        BigUint::from(rng.gen_range(1..curve.subgroup_order.to_u64_digits()[0]))
    } else {
        BigUint::from(args.known_log.unwrap_or(53))
    };
    if n >= curve.subgroup_order {
        return Err(format!(
            "known logarithm must be nonzero and smaller than subgroup order {}",
            curve.subgroup_order
        ));
    }
    Ok(n)
}
fn validate_factor_size(n: u32) -> Result<u32, String> {
    let dim = order_of_2_mod_n(n).ok_or("no supported factor-base family for this degree")?;
    if dim > MAX_FACTOR_DIMENSION {
        return Err(format!("factor-base dimension {dim} exceeds the materialization limit {MAX_FACTOR_DIMENSION}; inspection and fixture generation remain available"));
    }
    Ok(dim)
}
fn parameters(c: &KoblitzCurve, k: &BigUint, seed: u64) -> Parameters {
    let mut terms = c.curve.irreducible.low_terms.clone();
    terms.push(c.n);
    Parameters {
        schema_version: 1,
        name: format!("synthetic-k{}-n{}", c.a, c.n),
        field: Field::Binary {
            degree: c.n,
            polynomial_terms: terms,
        },
        a: c.a.to_string(),
        b: "1".into(),
        subgroup_order: c.subgroup_order.to_string(),
        cofactor: c.cofactor.to_string(),
        generator: params::binary_coordinates(c.generator()),
        point: params::binary_coordinates(&c.mul(c.generator(), k)),
        fixture: Some(Fixture {
            known_log: k.to_string(),
            seed,
        }),
    }
}
pub fn generate(args: GenerateArgs) -> Result<Value, String> {
    let c = curve(args.degree, args.curve_a)?;
    let run = RunArgs {
        degree: args.degree,
        curve_a: args.curve_a,
        random_target: true,
        seed: args.seed,
        ..RunArgs::default()
    };
    let k = known(&c, &run)?;
    serde_json::to_value(parameters(&c, &k, args.seed)).map_err(|e| e.to_string())
}
pub fn resources() -> Value {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        // libc defines the platform-specific layout; no shell/process inspection is used.
        let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) } == 0 {
            let multiplier = if cfg!(target_os = "macos") { 1 } else { 1024 };
            let seconds = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 / 1_000_000.0;
            return json!({"cpu_seconds":seconds(usage.ru_utime)+seconds(usage.ru_stime),
                "peak_rss_bytes":usage.ru_maxrss.max(0) as u64*multiplier,
                "scope":"process counters sampled before report emission; peak resident memory, not allocated bytes"});
        }
    }
    json!({"cpu_seconds":null,"peak_rss_bytes":null,"scope":"measurement unavailable on this platform"})
}
pub fn run(args: RunArgs, quiet: bool) -> Result<Value, String> {
    let begin = Instant::now();
    validate_factor_size(args.degree)?;
    if args.max_trials == 0 || args.max_trials > 20_000 {
        return Err("trial limit must be 1..=20000".into());
    }
    if !quiet {
        println!(
            "ic — synthetic index-calculus experiment\nConstructing K_{} / GF(2^{}) ...",
            args.curve_a, args.degree
        );
        let _ = std::io::stdout().flush();
    }
    let c = curve(args.degree, args.curve_a)?;
    let k = known(&c, &args)?;
    let supplied = parameters(&c, &k, args.seed);
    let target = c.mul(c.generator(), &k);
    if !quiet {
        println!(
            "Known logarithm: {k}; subgroup order: {}\nEngine: {}; factor candidate: {}",
            c.subgroup_order,
            args.solver.name(),
            args.factor_index
        );
    }
    let opts = KoblitzIcOptions {
        strategy: args.solver.strategy(),
        factor_index: args.factor_index,
        max_trials: args.max_trials as usize,
        seed: args.seed,
        collapse_negation: false,
        stop_on_verified_rank: false,
        allow_direct_relation: false,
        ..KoblitzIcOptions::default()
    };
    let mut stages = Vec::new();
    let mut stage_start = Instant::now();
    let report = koblitz_index_calculus_dlp_with_progress(&c, &target, &opts, &mut |event| {
        let (stage, state, details) = match event {
            KoblitzIcEvent::FactorBaseStarted => ("factor_base", "started", json!({})),
            KoblitzIcEvent::FactorBaseReady { points, orbits } => (
                "factor_base",
                "complete",
                json!({"points":points,"columns":orbits}),
            ),
            KoblitzIcEvent::RelationCollectionStarted { wanted } => {
                ("relation_collection", "started", json!({"wanted":wanted}))
            }
            KoblitzIcEvent::RelationProgress {
                collected,
                wanted,
                trials,
            } => {
                if !quiet {
                    println!(
                        "      {collected}/{wanted} relations; {trials} trials; {:.3}s",
                        stage_start.elapsed().as_secs_f64()
                    );
                    let _ = std::io::stdout().flush();
                }
                return;
            }
            KoblitzIcEvent::RelationCollectionFinished { collected, trials } => (
                "relation_collection",
                "stopped",
                json!({"collected":collected,"trials":trials}),
            ),
            KoblitzIcEvent::LinearAlgebraStarted { rows, columns } => (
                "linear_algebra",
                "started",
                json!({"rows":rows,"columns":columns}),
            ),
            KoblitzIcEvent::LinearAlgebraFinished => ("linear_algebra", "complete", json!({})),
            KoblitzIcEvent::LinearAlgebraIncomplete => ("linear_algebra", "incomplete", json!({})),
            KoblitzIcEvent::LinearAlgebraSkipped => ("linear_algebra", "skipped", json!({})),
            KoblitzIcEvent::VerificationStarted => ("verification", "started", json!({})),
            KoblitzIcEvent::VerificationFinished { verified } => (
                "verification",
                if verified { "pass" } else { "fail" },
                json!({}),
            ),
        };
        if state == "started" {
            stage_start = Instant::now();
        }
        let elapsed = stage_start.elapsed().as_secs_f64();
        if !quiet {
            let (index, label) = match stage {
                "factor_base" => (1, "Factor base"),
                "relation_collection" => (2, "Relation collection"),
                "linear_algebra" => (3, "Linear algebra"),
                _ => (4, "Verification"),
            };
            println!("[{index}/4] {label}: {state} {details} ({elapsed:.3}s)");
            let _ = std::io::stdout().flush();
        }
        stages.push(
            json!({"stage":stage,"status":state,"details":details,"elapsed_seconds":elapsed}),
        );
    });
    let Some(r) = report else {
        return Ok(
            json!({"schema_version":1,"operation":"run","status":"incomplete","evidence_scope":"synthetic_known_answer",
            "reason":"factor-base construction or pipeline operation did not complete; no result asserted",
            "arguments":args,"parameters":supplied,"stages":stages,"elapsed_seconds":begin.elapsed().as_secs_f64(),"resources":resources()}),
        );
    };
    let verified = r.sat_invalid_models == 0
        && r.log.as_ref() == Some(&k)
        && r.log
            .as_ref()
            .is_some_and(|d| c.mul(c.generator(), d) == target)
        && !r.direct_relation;
    Ok(
        json!({"schema_version":1,"operation":"run","status":if verified{"complete"}else{"incomplete"},
        "evidence_scope":"synthetic_known_answer","arguments":args,"parameters":supplied,"stages":stages,
        "result":{"expected":k.to_string(),"recovered":r.log.as_ref().map(ToString::to_string),"verified":verified},
        "counts":{"factor_base_points":r.factor_base_size,"columns":r.orbit_count,"relations":r.relations,"trials":r.trials,
            "f4_reductions":r.reductions,"sat_calls":r.sat_calls,"sat_unknowns":r.sat_unknowns,"sat_invalid_models":r.sat_invalid_models,
            "linear_solve_attempts":r.linear_solve_attempts,"cofactor_admissible":r.m_cofactor_admissible},
        "elapsed_seconds":begin.elapsed().as_secs_f64(),"resources":resources(),
        "limitations":["No imported target was used.","This run does not establish scaling or challenge readiness."]}),
    )
}
fn child(args: &RunArgs, seconds: u32) -> Result<Value, String> {
    let start = Instant::now();
    let mut command = Command::new(std::env::current_exe().map_err(|e| e.to_string())?);
    command.args([
        "run",
        "--json",
        "--degree",
        &args.degree.to_string(),
        "--curve-a",
        &args.curve_a.to_string(),
        "--factor-index",
        &args.factor_index.to_string(),
        "--seed",
        &args.seed.to_string(),
        "--max-trials",
        &args.max_trials.to_string(),
        "--solver",
        args.solver.name(),
        "--random-target",
    ]);
    let mut proc = command
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| e.to_string())?;
    let stdout = proc.stdout.take().unwrap();
    let stderr = proc.stderr.take().unwrap();
    let out = std::thread::spawn(move || {
        let mut v = Vec::new();
        stdout.take(1_048_577).read_to_end(&mut v).map(|_| v)
    });
    let err = std::thread::spawn(move || {
        let mut v = Vec::new();
        stderr.take(65_537).read_to_end(&mut v).map(|_| v)
    });
    let mut timed_out = false;
    let status = loop {
        if let Some(s) = proc.try_wait().map_err(|e| e.to_string())? {
            break s;
        }
        if start.elapsed() >= Duration::from_secs(seconds as u64) {
            timed_out = true;
            let _ = proc.kill();
            break proc.wait().map_err(|e| e.to_string())?;
        }
        std::thread::sleep(Duration::from_millis(10));
    };
    let stdout = out
        .join()
        .map_err(|_| "output reader failed")?
        .map_err(|e| e.to_string())?;
    let stderr = err
        .join()
        .map_err(|_| "error reader failed")?
        .map_err(|e| e.to_string())?;
    if timed_out {
        return Ok(
            json!({"status":"timed_out","elapsed_seconds":start.elapsed().as_secs_f64(),"arguments":args,"resources":null}),
        );
    }
    if stdout.len() > 1_048_576 {
        return Err("child report exceeded its output limit".into());
    }
    let mut report: Value = serde_json::from_slice(&stdout).unwrap_or_else(
        |_| json!({"status":"failed","error":String::from_utf8_lossy(&stderr),"arguments":args}),
    );
    if !status.success() && report["status"] == "complete" {
        report["status"] = json!("failed");
    }
    report["process_elapsed_seconds"] = json!(start.elapsed().as_secs_f64());
    Ok(report)
}
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 0 {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    }
}
pub fn compare(args: CompareArgs, quiet: bool) -> Result<Value, String> {
    let started = Instant::now();
    validate_factor_size(args.degree)?;
    let _ = curve(args.degree, args.curve_a)?;
    let count = factor_x_n_minus_1(args.degree).len();
    if count == 0 || count > 16 {
        return Err("candidate count is outside the bounded comparison range 1..=16".into());
    }
    let mut candidates = Vec::new();
    let mut winner = None;
    let mut best = f64::INFINITY;
    for index in 0..count {
        let mut trials = Vec::new();
        let mut costs = Vec::new();
        for sample in 0..args.samples {
            if !quiet {
                println!(
                    "Candidate {}/{}; training fixture {}/{}",
                    index + 1,
                    count,
                    sample + 1,
                    args.samples
                );
                let _ = std::io::stdout().flush();
            }
            let opts = RunArgs {
                degree: args.degree,
                curve_a: args.curve_a,
                random_target: true,
                seed: args.seed.wrapping_add(sample as u64),
                factor_index: index,
                max_trials: args.max_trials,
                solver: args.solver,
                ..RunArgs::default()
            };
            let result = child(&opts, args.timeout_seconds)?;
            if result["status"] == "complete" && result["result"]["verified"] == true {
                if let Some(cost) = result["process_elapsed_seconds"].as_f64() {
                    costs.push(cost);
                }
            }
            trials.push(result);
        }
        let eligible = costs.len() == args.samples as usize;
        let med = if eligible {
            Some(median(&mut costs))
        } else {
            None
        };
        if let Some(cost) = med {
            if cost < best {
                best = cost;
                winner = Some(index);
            }
        }
        candidates.push(json!({"factor_index":index,"eligible":eligible,"median_process_seconds":med,"training":trials}));
    }
    let mut holdout = Vec::new();
    if let Some(index) = winner {
        for sample in 0..args.holdout {
            if !quiet {
                println!(
                    "Candidate {index}; holdout fixture {}/{}",
                    sample + 1,
                    args.holdout
                );
                let _ = std::io::stdout().flush();
            }
            let opts = RunArgs {
                degree: args.degree,
                curve_a: args.curve_a,
                random_target: true,
                seed: (args.seed ^ 0x484f_4c44_4f55_5400).wrapping_add(sample as u64),
                factor_index: index,
                max_trials: args.max_trials,
                solver: args.solver,
                ..RunArgs::default()
            };
            holdout.push(child(&opts, args.timeout_seconds)?);
        }
    }
    let accepted = winner.is_some()
        && holdout.len() == args.holdout as usize
        && holdout
            .iter()
            .all(|v| v["status"] == "complete" && v["result"]["verified"] == true);
    Ok(
        json!({"schema_version":1,"operation":"compare","status":if accepted{"complete"}else{"inconclusive"},
        "evidence_scope":"bounded_synthetic_comparison","degree":args.degree,"curve_a":args.curve_a,
        "seed":args.seed,"solver":args.solver,"samples":args.samples,"holdout_samples":args.holdout,
        "candidate_family":"degree-ord_n(2) irreducible factors used by the existing materialized builder",
        "candidate_count":count,"candidates":candidates,"training_winner":winner,"selected_factor_index":if accepted{winner}else{None},
        "holdout":holdout,"elapsed_seconds":started.elapsed().as_secs_f64(),"resources":resources(),
        "selection_metric":"median child-process wall time including startup, construction, collection, solve, verification, and report emission",
        "completion_poll_interval_ms":10,
        "scope":"fastest observed eligible candidate on these training fixtures, with correctness checked on separate holdout fixtures",
        "limitations":["Not a global optimum or an asymptotic result.","Failed and timed-out cases remain in the report and are ineligible.",
            "Child resource metrics cover separate processes; parent metrics exclude child CPU and memory.","No imported target was used."]}),
    )
}
