//! `ic bench` — run index-calculus configurations and compare them.
//!
//! ```text
//! ic bench --list                       what can be plugged in where
//! ic bench --bits 16 --factor-base prime-abscissa:size=24 \
//!          --oracle mitm:negation_folded=1 --targets walk
//! ic bench --sweep sweep.json --out docs/ic/runs/sweep-YYYY-MM-DD.json
//! ```
//!
//! A configuration is a choice at each stage.  The table it prints
//! carries every stage's numbers, so the effect of a choice is visible
//! where it lands: change the factor base and the hit rate moves, which
//! moves the trials, which moves the matrix.
//!
//! See `docs/ic/FRAMEWORK.md`.

use clap::Args;
use serde_json::{json, Value};

use crypto_lib::cryptanalysis::ic_boundary::{
    calibrate_group, calibrate_row_ops, koblitz_instance, random_binary_instance, roster_prime_instance, BinaryGroup,
    BinaryInstance, Calibration, CountedGroup, GroupOps, PrimeInstance,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::{
    BinarySubspaceBase, FrobeniusMitmOracle, KoblitzOrbitBase, MitmOracle, PrimeAbscissaBase,
    SubtractOracle,
};
use crypto_lib::cryptanalysis::ic_framework::solvers::{solver_by_name, solver_registry};
use crypto_lib::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, Targets,
};
use crypto_lib::cryptanalysis::ic_framework::{format_markdown, run_pipeline, PipelineSpec, RunReport};

#[derive(Args, Clone)]
pub struct BenchArgs {
    /// Print what can be plugged in at each stage, with the parameters
    /// each plug-in reads, and exit.
    #[arg(long)]
    pub list: bool,
    /// Prime-field instance of this subgroup bit length.
    #[arg(long)]
    pub bits: Option<u32>,
    /// Binary-field instance of this field degree (random curve).
    #[arg(long)]
    pub char2_degree: Option<u32>,
    /// Koblitz instance of this field degree.
    #[arg(long)]
    pub koblitz_degree: Option<u32>,
    /// `name` or `name:k=v,k=v`.
    #[arg(long)]
    pub factor_base: Option<String>,
    /// `name` or `name:k=v,k=v`.
    #[arg(long)]
    pub oracle: Option<String>,
    /// The polynomial solver, for oracles that use one.
    #[arg(long)]
    pub solver: Option<String>,
    /// `random` or `walk`.
    #[arg(long, default_value = "walk")]
    pub targets: String,
    #[arg(long, default_value_t = 2_000_000)]
    pub max_trials: u64,
    #[arg(long, default_value_t = 0x1C_B0_0B_DA_7A)]
    pub seed: u64,
    /// Repeats per configuration; each draws a fresh target.
    #[arg(long, default_value_t = 1)]
    pub repeats: usize,
    /// Per-solver-call budget in seconds; 0 for none.
    #[arg(long, default_value_t = 120)]
    pub solver_budget_seconds: u64,
    /// A JSON file of configurations to sweep.  See the framework
    /// documentation for the schema.
    #[arg(long)]
    pub sweep: Option<std::path::PathBuf>,
}

/// `name:k=v,k=v` → `(name, params)`.
fn parse_plugin(spec: &str) -> Result<(String, Params), String> {
    let (name, rest) = match spec.split_once(':') {
        None => (spec, ""),
        Some((n, r)) => (n, r),
    };
    let mut params = Params::default();
    for kv in rest.split(',').filter(|s| !s.is_empty()) {
        let (k, v) = kv
            .split_once('=')
            .ok_or_else(|| format!("parameter `{kv}` is not key=value"))?;
        params.set(k.trim(), v.trim());
    }
    Ok((name.to_string(), params))
}

fn listing() -> Value {
    let solvers: Vec<Value> = solver_registry()
        .iter()
        .map(|s| {
            json!({
                "name": s.name(),
                "describes": s.describe(),
                "parameters": s.parameters().iter()
                    .map(|(k, d)| json!({"name": k, "means": d}))
                    .collect::<Vec<_>>(),
            })
        })
        .collect();
    json!({
        "stages": [
            {
                "stage": "factor base",
                "trait": "FactorBaseBuilder",
                "what_it_chooses": "which points relations are written over, and which unknown each point contributes to",
                "plugins": [
                    {"name": "prime-abscissa", "regimes": ["prime"],
                     "parameters": [{"name": "size", "means": "abscissae in the base; the family optimum is about #E^(1/3)"}]},
                    {"name": "binary-subspace", "regimes": ["char2", "koblitz"],
                     "parameters": [{"name": "dimension", "means": "F_2-dimension of the abscissa subspace"}]},
                    {"name": "koblitz-orbit", "regimes": ["koblitz"],
                     "parameters": [
                        {"name": "divisor", "means": "comma-separated factor indices selecting the Frobenius-invariant subspace"},
                        {"name": "no_fold", "means": "1 for one column per abscissa: the control that shows what the fold buys"}]},
                ],
            },
            {
                "stage": "targets",
                "trait": "Targets",
                "what_it_chooses": "how each trial point is produced",
                "plugins": [
                    {"name": "random", "means": "a fresh [a]G + [b]Q per trial: two scalar multiplications"},
                    {"name": "walk", "means": "an r-adding walk: one addition per trial, with a repeat guard"},
                ],
            },
            {
                "stage": "point decomposition",
                "trait": "DecompositionOracle",
                "what_it_chooses": "how a target is written as a sum of factor-base points",
                "plugins": [
                    {"name": "subtract", "summands": 2, "regimes": ["prime", "char2", "koblitz"],
                     "parameters": []},
                    {"name": "mitm", "summands": "2 or 3", "regimes": ["prime", "char2", "koblitz"],
                     "parameters": [{"name": "negation_folded", "means": "1 to halve the table's additions"},
                                    {"name": "m", "means": "summands, 2 or 3"}]},
                    {"name": "mitm-frobenius", "summands": "2 or 3", "regimes": ["koblitz"],
                     "parameters": [{"name": "m", "means": "summands, 2 or 3"}]},
                ],
            },
            {
                "stage": "polynomial system solver",
                "trait": "SystemSolver",
                "what_it_chooses": "how an algebraic oracle decides its system; the plug point for F4, F5, XL, SAT",
                "plugins": solvers,
            },
            {
                "stage": "relation matrix",
                "trait": "RelationSolver",
                "what_it_chooses": "how relations are accumulated and the logarithm read off",
                "plugins": [
                    {"name": "incremental-gauss",
                     "means": "reduced row echelon maintained as rows arrive, stopping the moment the target column is pinned"},
                ],
            },
        ],
        "how_to_add_one": "implement the stage's trait and register it; docs/ic/FRAMEWORK.md walks through a solver end to end",
    })
}

/// Measure this host, then pin every ratio the repository's table
/// carries, so two runs of the same counts price the same (§12 of the
/// ledger note).
fn calibration_for<G: CountedGroup>(
    g: &G,
    points: &[G::Elt],
    regime: &str,
    instance: &str,
    modulus: u64,
) -> Calibration {
    let mut calib = Calibration::default();
    if !points.is_empty() {
        calibrate_group(g, points, &mut calib);
    }
    calibrate_row_ops(modulus, &mut calib);
    calib.pin(regime, instance);
    calib
}

fn sample_prime_points(inst: &PrimeInstance) -> Vec<crypto_lib::cryptanalysis::ic_boundary::PrimePoint> {
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    (1..=8u64).map(|k| inst.curve.mul(&mut ops, g, k)).collect()
}

fn sample_binary_points(
    inst: &BinaryInstance,
) -> Vec<crypto_lib::cryptanalysis::koblitz_fast::FastPoint> {
    let g = BinaryGroup(&inst.fast);
    let mut ops = GroupOps::default();
    (1..=8u64).map(|k| g.mul(&mut ops, inst.generator, k)).collect()
}

/// One prime-field configuration, over `repeats` targets.
fn run_prime(
    inst: &PrimeInstance,
    spec: &PipelineSpec,
    args: &BenchArgs,
    calib: &Calibration,
) -> Result<Vec<RunReport>, String> {
    let (fb_name, _) = parse_plugin(args.factor_base.as_deref().unwrap_or("prime-abscissa"))?;
    let (or_name, or_params) = parse_plugin(args.oracle.as_deref().unwrap_or("mitm"))?;
    if fb_name != "prime-abscissa" {
        return Err(format!(
            "factor base `{fb_name}` is not available on a prime-field curve; try prime-abscissa"
        ));
    }
    let base = PrimeAbscissaBase { instance: inst };
    let mut out = Vec::new();
    for rep in 0..args.repeats.max(1) {
        let mut ops = GroupOps::default();
        let g = inst.generator_point();
        let planted = 1 + (spec.seed.wrapping_add(rep as u64 * 0x9E37)) % (inst.r - 1);
        let q = inst.curve.mul(&mut ops, g, planted);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: g,
            target: q,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let mut spec = spec.clone();
        spec.seed = spec.seed.wrapping_add(rep as u64 * 0x9E37);
        let m = or_params.u64_or("m", 2)? as u32;
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(m);
        let oracle: &mut dyn DecompositionOracle<_> = match or_name.as_str() {
            "subtract" => &mut subtract,
            "mitm" => &mut mitm,
            other => {
                return Err(format!(
                    "oracle `{other}` is not available on a prime-field curve; try subtract or mitm"
                ))
            }
        };
        out.push(run_pipeline(
            &ctx, &spec, &base, oracle, None, planted, calib, None,
        )?);
    }
    Ok(out)
}

/// One binary or Koblitz configuration, over `repeats` targets.
fn run_binary(
    inst: &BinaryInstance,
    spec: &PipelineSpec,
    args: &BenchArgs,
    calib: &Calibration,
) -> Result<Vec<RunReport>, String> {
    let (fb_name, _) = parse_plugin(args.factor_base.as_deref().unwrap_or("binary-subspace"))?;
    let (or_name, or_params) = parse_plugin(args.oracle.as_deref().unwrap_or("mitm"))?;
    let g = BinaryGroup(&inst.fast);
    let subspace = BinarySubspaceBase { instance: inst };
    let orbit = KoblitzOrbitBase { instance: inst };
    let base: &dyn FactorBaseBuilder<BinaryGroup> = match fb_name.as_str() {
        "binary-subspace" => &subspace,
        "koblitz-orbit" => &orbit,
        other => {
            return Err(format!(
                "factor base `{other}` is not available on a binary curve; try binary-subspace or koblitz-orbit"
            ))
        }
    };
    let mut out = Vec::new();
    for rep in 0..args.repeats.max(1) {
        let mut ops = GroupOps::default();
        let planted = 1 + (spec.seed.wrapping_add(rep as u64 * 0x9E37)) % (inst.r - 1);
        let q = g.mul(&mut ops, inst.generator, planted);
        let ctx = InstanceCtx {
            group: &g,
            generator: inst.generator,
            target: q,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        let mut spec = spec.clone();
        spec.seed = spec.seed.wrapping_add(rep as u64 * 0x9E37);
        let m = or_params.u64_or("m", 2)? as u32;
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(m);
        let mut frob = FrobeniusMitmOracle::new(m, inst);
        let oracle: &mut dyn DecompositionOracle<BinaryGroup> = match or_name.as_str() {
            "subtract" => &mut subtract,
            "mitm" => &mut mitm,
            "mitm-frobenius" => &mut frob,
            other => return Err(format!("oracle `{other}` is not available here")),
        };
        out.push(run_pipeline(
            &ctx, &spec, base, oracle, None, planted, calib, None,
        )?);
    }
    Ok(out)
}

/// The instance a sweep runs every configuration on, built once.
enum Instance {
    Prime(PrimeInstance),
    Binary(BinaryInstance),
}

impl Instance {
    fn name(&self) -> &str {
        match self {
            Instance::Prime(i) => &i.name,
            Instance::Binary(i) => &i.name,
        }
    }
}

/// A sweep file: one instance, and either an explicit list of
/// configurations or a matrix to take the product of.
#[derive(serde::Deserialize)]
struct SweepFile {
    instance: SweepInstance,
    #[serde(default)]
    repeats: Option<usize>,
    #[serde(default)]
    seed: Option<u64>,
    #[serde(default)]
    max_trials: Option<u64>,
    /// Explicit configurations, each a map of stage to plug-in spec.
    #[serde(default)]
    configurations: Vec<std::collections::BTreeMap<String, String>>,
    /// A matrix: every combination of the listed choices is run.
    #[serde(default)]
    matrix: std::collections::BTreeMap<String, Vec<String>>,
}

#[derive(serde::Deserialize)]
struct SweepInstance {
    regime: String,
    degree: u32,
}

/// Every combination of the matrix, in a stable order.
fn expand_matrix(
    matrix: &std::collections::BTreeMap<String, Vec<String>>,
) -> Vec<std::collections::BTreeMap<String, String>> {
    let mut out = vec![std::collections::BTreeMap::new()];
    for (key, choices) in matrix {
        if choices.is_empty() {
            continue;
        }
        let mut next = Vec::with_capacity(out.len() * choices.len());
        for base in &out {
            for choice in choices {
                let mut row = base.clone();
                row.insert(key.clone(), choice.clone());
                next.push(row);
            }
        }
        out = next;
    }
    out
}

fn spec_from(
    cfg: &std::collections::BTreeMap<String, String>,
    defaults: &BenchArgs,
) -> Result<(PipelineSpec, String, String), String> {
    let fb = cfg
        .get("factor_base")
        .cloned()
        .or_else(|| defaults.factor_base.clone())
        .ok_or("a configuration needs a `factor_base`")?;
    let or = cfg
        .get("oracle")
        .cloned()
        .or_else(|| defaults.oracle.clone())
        .ok_or("a configuration needs an `oracle`")?;
    let tg = cfg
        .get("targets")
        .cloned()
        .unwrap_or_else(|| defaults.targets.clone());
    let (fb_name, fb_params) = parse_plugin(&fb)?;
    let (or_name, or_params) = parse_plugin(&or)?;
    Ok((
        PipelineSpec {
            factor_base: fb_name,
            factor_base_params: fb_params,
            oracle: or_name,
            oracle_params: or_params,
            solver: cfg.get("solver").cloned().or_else(|| defaults.solver.clone()),
            solver_params: Params::default(),
            targets: Targets::parse(&tg)?,
            max_trials: defaults.max_trials,
            seed: defaults.seed,
            solver_budget_seconds: defaults.solver_budget_seconds,
        },
        fb,
        or,
    ))
}

pub fn run(args: BenchArgs, json_only: bool) -> Result<Value, String> {
    if args.list {
        let l = listing();
        if !json_only {
            eprintln!("{}", serde_json::to_string_pretty(&l).unwrap_or_default());
        }
        return Ok(json!({
            "schema_version": 1,
            "operation": "bench-list",
            "status": "complete",
            "registry": l,
        }));
    }

    // Validate the solver name up front, so a typo fails before a long
    // run rather than after it.
    if let Some(name) = &args.solver {
        let (n, _) = parse_plugin(name)?;
        solver_by_name(&n)?;
    }

    // Either a sweep file or the flags describe the work.
    let (instance_spec, configs, mut args) = match &args.sweep {
        Some(path) => {
            let raw = std::fs::read_to_string(path)
                .map_err(|e| format!("cannot read sweep file {}: {e}", path.display()))?;
            let file: SweepFile = serde_json::from_str(&raw)
                .map_err(|e| format!("sweep file {} is not valid: {e}", path.display()))?;
            let mut a = args.clone();
            if let Some(v) = file.repeats {
                a.repeats = v;
            }
            if let Some(v) = file.seed {
                a.seed = v;
            }
            if let Some(v) = file.max_trials {
                a.max_trials = v;
            }
            let mut configs = file.configurations.clone();
            configs.extend(expand_matrix(&file.matrix));
            if configs.is_empty() {
                return Err("the sweep file lists no configurations and no matrix".into());
            }
            ((file.instance.regime.clone(), file.instance.degree), configs, a)
        }
        None => {
            let (regime, degree) = match (args.bits, args.char2_degree, args.koblitz_degree) {
                (Some(b), None, None) => ("prime".to_string(), b),
                (None, Some(n), None) => ("char2".to_string(), n),
                (None, None, Some(n)) => ("koblitz".to_string(), n),
                _ => {
                    return Err(
                        "choose exactly one of --bits, --char2-degree, --koblitz-degree, or --sweep"
                            .into(),
                    )
                }
            };
            let mut one = std::collections::BTreeMap::new();
            if let Some(v) = &args.factor_base {
                one.insert("factor_base".to_string(), v.clone());
            }
            if let Some(v) = &args.oracle {
                one.insert("oracle".to_string(), v.clone());
            }
            one.insert("targets".to_string(), args.targets.clone());
            ((regime, degree), vec![one], args.clone())
        }
    };

    let (regime, degree) = instance_spec;
    let instance = match regime.as_str() {
        "prime" => Instance::Prime(
            roster_prime_instance(degree)
                .ok_or_else(|| format!("no prime instance at {degree} bits"))?,
        ),
        "char2" => Instance::Binary(
            random_binary_instance(degree, args.seed, 1 << 20)
                .ok_or_else(|| format!("no random binary instance at degree {degree}"))?,
        ),
        "koblitz" => Instance::Binary(
            koblitz_instance(1, degree)
                .or_else(|| koblitz_instance(0, degree))
                .ok_or_else(|| format!("no Koblitz instance at degree {degree}"))?,
        ),
        other => return Err(format!("unknown regime `{other}`; try prime, char2 or koblitz")),
    };

    let calib = match &instance {
        Instance::Prime(i) => {
            calibration_for(&i.curve, &sample_prime_points(i), "prime", &i.name, i.r)
        }
        Instance::Binary(i) => {
            let g = BinaryGroup(&i.fast);
            calibration_for(&g, &sample_binary_points(i), &regime, &i.name, i.r)
        }
    };

    let started = std::time::Instant::now();
    let mut rows: Vec<RunReport> = Vec::new();
    let mut failures: Vec<Value> = Vec::new();
    for cfg in &configs {
        let (spec, fb, or) = spec_from(cfg, &args)?;
        args.factor_base = Some(fb);
        args.oracle = Some(or);
        if !json_only {
            eprintln!("  {} …", spec.label());
        }
        // A configuration that does not apply to this regime is
        // reported, not fatal: a sweep over a matrix will contain
        // combinations that do not exist, and the useful output is the
        // ones that do plus a note on the ones that do not.
        let outcome = match &instance {
            Instance::Prime(i) => run_prime(i, &spec, &args, &calib),
            Instance::Binary(i) => run_binary(i, &spec, &args, &calib),
        };
        match outcome {
            Ok(mut got) => rows.append(&mut got),
            Err(why) => failures.push(json!({"configuration": spec.label(), "why": why})),
        }
    }

    let all_verified = !rows.is_empty() && rows.iter().all(|r| r.verified);
    let markdown = format_markdown(&rows);
    if !json_only {
        eprintln!("\n{markdown}");
        for f in &failures {
            eprintln!("  skipped: {}", f["configuration"].as_str().unwrap_or("?"));
        }
    }
    Ok(json!({
        "schema_version": 1,
        "operation": "bench",
        "status": if all_verified { "complete" } else { "incomplete" },
        "what_this_is": "Index-calculus configurations -- each a choice of factor base, target source, decomposition oracle, polynomial solver and relation matrix -- run end to end against a planted logarithm, with every stage's cost in group-addition equivalents and the total divided by sqrt(r).",
        "what_this_is_not": [
            "not a speed claim: a configuration on a toy instance is a measurement, and a comparison needs the matched reference rows AGENTS.md section 1 asks for",
            "not a wall-clock benchmark: operation counts are the metric, wall time rides along",
            "not a claim about any deployed curve"
        ],
        "instance": instance.name(),
        "regime": regime,
        "configurations_run": rows.len(),
        "configurations_skipped": failures,
        "all_verified": all_verified,
        "elapsed_seconds": started.elapsed().as_secs_f64(),
        "markdown": markdown,
        "rows": rows,
    }))
}
