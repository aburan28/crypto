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

use crypto_lib::cryptanalysis::glv_invariant_base::{generate_cm_instance, CmFamily};
use crypto_lib::cryptanalysis::ic_boundary::{
    calibrate_binary_instance, calibrate_group, calibrate_row_ops, calibrate_word_xor,
    generic_floor_ops, koblitz_instance, random_binary_instance, rho_reference,
    rho_reference_negation, roster_prime_instance, signed_frobenius_rho, BinaryGroup,
    BinaryInstance, Calibration, CountedGroup, GroupOps, PinOutcome, PrimeInstance, RhoResult,
};
use crypto_lib::cryptanalysis::ic_framework::linalg::MATRIX_NAMES;
use crypto_lib::cryptanalysis::ic_framework::plugins::{
    BinarySubspaceBase, DescentAlgebraicOracle, FrobeniusMitmOracle, GlvOrbitBase,
    KoblitzOrbitBase, KoblitzSymmetrisedBase, KoblitzTraceZeroBase, MitmOracle, PrimeAbscissaBase,
    SubtractOracle, SymmetrisedOracle,
};
use crypto_lib::cryptanalysis::ic_framework::solvers::{
    solver_by_name, solver_registry, validate_solver_params,
};
use crypto_lib::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, Targets,
};
use crypto_lib::cryptanalysis::ic_framework::{
    format_markdown, run_pipeline, PipelineSpec, RunReport,
};

#[derive(Args, Clone)]
pub struct BenchArgs {
    /// Print what can be plugged in at each stage, with the parameters
    /// each plug-in reads, and exit.
    #[arg(long)]
    pub list: bool,
    /// Prime-field instance of this subgroup bit length.
    #[arg(long)]
    pub bits: Option<u32>,
    /// With `--bits`: the prime-field curve family.  `generic` takes the
    /// roster's curve at that size; `j0`, `j1728`, `d7` and `d8` generate
    /// a CM curve of that family (order certified against the Cornacchia
    /// candidates, deterministic in `--seed`) so that `glv-orbit` has an
    /// automorphism group to fold by.
    #[arg(long, default_value = "generic")]
    pub family: String,
    /// Binary-field instance of this field degree (random curve).
    #[arg(long)]
    pub char2_degree: Option<u32>,
    /// Koblitz instance of this field degree.
    #[arg(long)]
    pub koblitz_degree: Option<u32>,
    /// Which Koblitz curve at that degree: `0` for `K_0` (`a = 0`) or `1`
    /// for `K_1`.  Without it the bench takes `K_1` and falls back to
    /// `K_0`, so a ladder over both curves at one degree needs this.
    #[arg(long)]
    pub koblitz_a: Option<u8>,
    /// Largest cofactor a random binary curve may have.  The boundary
    /// ladders use 8, so `S = ops / sqrt(r)` is taken over a subgroup
    /// close to the whole group, as on every ledger row; a larger bound
    /// finds a curve faster but can hand back a tiny `r` that makes `S`
    /// meaningless against rho.
    #[arg(long, default_value_t = 8)]
    pub max_cofactor: u64,
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
    /// The relation matrix: `incremental-gauss` or `structured-gauss`.
    #[arg(long, default_value = "incremental-gauss")]
    pub linalg: String,
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
    /// Counted Pollard-rho runs on the same instance and planted targets,
    /// averaged into the `vs rho` column; 0 leaves the column empty.
    #[arg(long, default_value_t = 16)]
    pub rho_runs: usize,
    /// A JSON file of configurations to sweep.  See the framework
    /// documentation for the schema.
    #[arg(long)]
    pub sweep: Option<std::path::PathBuf>,
    /// Write every factor base the run builds, point by point, to this
    /// JSON file (a new file; an existing one is never overwritten), so
    /// the base a report's rows were taken over can be frozen next to
    /// the report.  Binary and Koblitz instances only.
    #[arg(long)]
    pub factor_base_out: Option<std::path::PathBuf>,
}

/// The factor base a configuration builds, frozen point by point: the
/// curve it was built on, the builder and its parameters, the column
/// each point folds onto, and a hash of the point list so two reports
/// can be checked against the same base.
fn frozen_factor_base<'c>(
    inst: &BinaryInstance,
    spec: &PipelineSpec,
    base: &dyn FactorBaseBuilder<BinaryGroup<'c>>,
    ctx: &InstanceCtx<'_, BinaryGroup<'c>>,
) -> Result<Value, String> {
    let mut ops = GroupOps::default();
    let fb = base.build(ctx, &spec.factor_base_params, &mut ops)?;
    let f = &inst.fast.field;
    let mut hasher = blake3::Hasher::new();
    let points: Vec<Value> = fb
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            hasher.update(&p.x.to_le_bytes());
            hasher.update(&p.y.to_le_bytes());
            // The frame the symmetrised system works in: u = 1/(x + 1),
            // in which the rational 2-torsion acts as u ↦ u + 1.
            let u = (inst.koblitz.is_some() && p.x != 1).then(|| f.inv(p.x ^ 1));
            json!({
                "i": i,
                "x": format!("{:#x}", p.x),
                "y": format!("{:#x}", p.y),
                "u": u.map(|u| format!("{u:#x}")),
                "col": fb.col_of.get(i),
                "coef": fb.coef_of.get(i),
                "neg": fb.neg_index.get(i),
            })
        })
        .collect();
    Ok(json!({
        "schema_version": 1,
        "what_this_is": "The factor base one `ic bench` configuration was run over, frozen point by point on the curve it was built on.",
        "instance": inst.name,
        "curve": {
            "degree": inst.n,
            "irreducible_low_terms": inst.irreducible.low_terms,
            "a": format!("{:#x}", inst.a),
            "b": format!("{:#x}", inst.b),
            "koblitz_a": inst.koblitz.as_ref().map(|k| k.a),
            "group_order": inst.group_order,
            "r": inst.r,
            "cofactor": inst.cofactor,
            "generator": {"x": format!("{:#x}", inst.generator.x), "y": format!("{:#x}", inst.generator.y)},
        },
        "factor_base": spec.factor_base,
        "factor_base_params": spec.factor_base_params,
        "description": fb.description,
        "dimension": fb.dimension,
        "columns": fb.columns,
        "abscissae": fb.abscissae,
        "points": fb.points.len(),
        "subspace_basis": fb.subspace_basis.as_ref().map(|b| b.iter().map(|w| format!("{w:#x}")).collect::<Vec<_>>()),
        "build_cost": fb.cost,
        "point_list_blake3": hasher.finalize().to_hex().to_string(),
        "point_list": points,
    }))
}

/// `name:k=v,k=v` → `(name, params)`; list values use `;`.
fn parse_plugin(spec: &str) -> Result<(String, Params), String> {
    Params::parse_spec(spec)
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
                    {"name": "glv-orbit", "regimes": ["prime"],
                     "parameters": [
                        {"name": "size", "means": "seed abscissae; the base is their closure under the curve's automorphism group"},
                        {"name": "group", "means": "auto, negation, j0 or j1728 (use --family j0 or j1728 for a curve that has one)"},
                        {"name": "no_fold", "means": "1 to fold the same points by negation only: the control that shows what the fold buys"}]},
                    {"name": "gls-line", "regimes": ["gls (examples/glv_invariant_bench.rs)"],
                     "parameters": [
                        {"name": "no_fold", "means": "1 to fold the same points by negation only"}]},
                    {"name": "binary-subspace", "regimes": ["char2", "koblitz"],
                     "parameters": [{"name": "dimension", "means": "F_2-dimension of the abscissa subspace"}]},
                    {"name": "koblitz-orbit", "regimes": ["koblitz"],
                     "parameters": [
                        {"name": "divisor", "means": "`;`-separated factor indices selecting the Frobenius-invariant subspace, e.g. divisor=1;2 (a `,` would end the parameter)"},
                        {"name": "no_fold", "means": "1 for one column per abscissa: the control that shows what the fold buys"}]},
                    {"name": "koblitz-trace-zero", "regimes": ["koblitz"],
                     "parameters": [
                        {"name": "divisor", "means": "`;`-separated factor indices; the selected invariant abscissa subspace must lie entirely in the absolute-trace kernel"}]},
                    {"name": "koblitz-symmetrised", "regimes": ["koblitz"],
                     "parameters": [
                        {"name": "divisor", "means": "`;`-separated factor indices, e.g. divisor=0;1; must include 0 (x + 1) so that 1 ∈ V; the base is F_u = {P : 1/(x+1) ∈ V}, T = (0,1) included"},
                        {"name": "no_fold", "means": "1 for one column per abscissa"}]},
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
                    {"name": "descent-algebraic", "summands": "2 or 3", "regimes": ["char2", "koblitz"],
                     "parameters": [{"name": "m", "means": "summands, 2 (descends S3) or 3 (descends S4)"}],
                     "needs": "a subspace factor base (binary-subspace or koblitz-orbit) and --solver"},
                    {"name": "symmetrised", "summands": "2 or 3", "regimes": ["koblitz"],
                     "parameters": [{"name": "m", "means": "summands, 2 (symmetrised S3) or 3 (symmetrised S4)"},
                                    {"name": "divisor", "means": "the base's divisor; copied from koblitz-symmetrised when omitted"},
                                    {"name": "engine", "means": "inherited-f4 (default), matrix-f4 or matrix-f5"},
                                    {"name": "max_degree", "means": "Macaulay cap before splitting (default 3)"},
                                    {"name": "node_budget", "means": "algebraic reduction calls before a solve gives up (default 4096)"}],
                     "needs": "the koblitz-symmetrised factor base; solves in w = u² + u, s = Σu with koblitz_groebner"},
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
                "plugins": MATRIX_NAMES.iter().map(|(n, d)| json!({"name": n, "means": d})).collect::<Vec<_>>(),
            },
        ],
        "how_to_add_one": "implement the stage's trait and register it; docs/ic/FRAMEWORK.md walks through a solver end to end",
    })
}

/// Measure this host, then pin every ratio the repository's table
/// carries, so two runs of the same counts price the same (§12 of the
/// ledger note).
/// The prime regime's calibration: the group, the matrix and the word
/// XOR measured here, the prime-only units (square roots, Legendre
/// symbols, inversions) taken from the pinned table for the roster
/// curves, and every unit the table carries pinned over the measured
/// value.
fn calibration_for<G: CountedGroup>(
    g: &G,
    points: &[G::Elt],
    regime: &str,
    instance: &str,
    modulus: u64,
) -> (Calibration, PinOutcome) {
    let mut calib = Calibration::default();
    if !points.is_empty() {
        calibrate_group(g, points, &mut calib);
    }
    calibrate_row_ops(modulus, &mut calib);
    calib.ns_per_word_xor = Some(calibrate_word_xor());
    let pins = calib.pin(regime, instance);
    (calib, pins)
}

/// The binary regimes' calibration: the ledger's own measurement of
/// the group, the Artin–Schreier solve, the Frobenius map, the matrix
/// and the word XOR, then the pinned ratios where the table has this
/// instance.  A freshly generated random curve is not in the table, so
/// its base build is priced at the measured Artin–Schreier factor
/// rather than left at zero.
fn calibration_for_binary(inst: &BinaryInstance, regime: &str) -> (Calibration, PinOutcome) {
    let mut calib = calibrate_binary_instance(inst);
    calibrate_row_ops(inst.r, &mut calib);
    let pins = calib.pin(regime, &inst.name);
    (calib, pins)
}

fn sample_prime_points(
    inst: &PrimeInstance,
) -> Vec<crypto_lib::cryptanalysis::ic_boundary::PrimePoint> {
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    (1..=8u64).map(|k| inst.curve.mul(&mut ops, g, k)).collect()
}

/// The logarithm a configuration's `rep`-th repeat plants.
fn planted_log(seed: u64, rep: usize, r: u64) -> u64 {
    1 + (seed.wrapping_add(rep as u64 * 0x9E37)) % (r - 1)
}

/// One counted rho run.
#[derive(serde::Serialize)]
struct RhoRun {
    planted: u64,
    seed: u64,
    s: f64,
    /// Walk operations over `√r`: the part the `√(πr/2A)` floor is about.
    s_walk: f64,
    steps: u64,
    steps_over_expected: f64,
    verified: bool,
}

/// **The method boundary, measured on the instance.**  A counted Pollard
/// rho run on the logarithms the configurations plant, several walks
/// averaged because one rho run's cost is a draw from a wide
/// distribution.
#[derive(serde::Serialize)]
struct RhoReference {
    method: String,
    /// Automorphisms the walk uses: `2` for the negation map, `2n` for
    /// the signed Frobenius, `1` for the plain walk.
    automorphisms: u32,
    runs: usize,
    mean_s: f64,
    min_s: f64,
    max_s: f64,
    mean_s_walk: f64,
    all_verified: bool,
    per_run: Vec<RhoRun>,
}

/// `runs` rho runs of one walk on the configurations' planted
/// logarithms.  `walk(target, planted, seed, max_steps)` runs one; the
/// seeds and targets are the same for every walk, so two references on
/// one instance are paired.
fn rho_reference_for<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    r: u64,
    seed: u64,
    repeats: usize,
    runs: usize,
    mut walk: impl FnMut(G::Elt, u64, u64, u64) -> RhoResult,
) -> Option<RhoReference> {
    if runs == 0 || r < 3 {
        return None;
    }
    let max_steps = (generic_floor_ops(r as f64, 1.0) * 64.0) as u64 + 4096;
    let mut per_run = Vec::with_capacity(runs);
    let (mut method, mut automorphisms) = (String::new(), 1);
    for k in 0..runs {
        let planted = planted_log(seed, k % repeats.max(1), r);
        let mut ops = GroupOps::default();
        let target = g.mul(&mut ops, generator, planted);
        let run_seed = seed ^ (0x5248_4F00 + k as u64);
        let res = walk(target, planted, run_seed, max_steps);
        method = res.method.clone();
        automorphisms = res.automorphisms;
        per_run.push(RhoRun {
            planted,
            seed: run_seed,
            s: res.s,
            s_walk: res.s_walk,
            steps: res.steps,
            steps_over_expected: res.steps_over_expected,
            verified: res.verified && res.recovered == Some(planted),
        });
    }
    let s: Vec<f64> = per_run.iter().map(|r| r.s).collect();
    Some(RhoReference {
        method,
        automorphisms,
        runs,
        mean_s: s.iter().sum::<f64>() / runs as f64,
        min_s: s.iter().copied().fold(f64::INFINITY, f64::min),
        max_s: s.iter().copied().fold(0.0, f64::max),
        mean_s_walk: per_run.iter().map(|r| r.s_walk).sum::<f64>() / runs as f64,
        all_verified: per_run.iter().all(|r| r.verified),
        per_run,
    })
}

/// The references an instance is priced against.
struct RhoReferences {
    /// The matched reference: the `vs rho` column divides by its mean.
    matched: Option<RhoReference>,
    /// The plain walk (`A = 1`) every run was priced against through
    /// ledger §17, on the same seeds: the before mark.
    plain: Option<RhoReference>,
    /// The other eligible automorphism-aware walks that were run on the
    /// same seeds (Koblitz curves), kept beside the one chosen.
    candidates: Vec<RhoReference>,
    /// How the matched reference was chosen.
    rule: &'static str,
}

/// **The matched rho** (accounting contract, `comparison_contract.rho`):
/// the negation map on the prime and random binary curves; on a Koblitz
/// curve both eligible walks — the repository's signed-Frobenius walk
/// (`A = 2n`) and the negation walk — and the one with the lower mean
/// `S`, because on toy subgroups the signed walk's set-up (thirty-two
/// parallel walks started by scalar multiplications) outweighs what
/// the Frobenius saves.  The plain walk rides along as the before mark.
fn matched_rho(instance: &Instance, seed: u64, repeats: usize, runs: usize) -> RhoReferences {
    match instance {
        Instance::Prime(i) => {
            let (g, gen, r) = (&i.curve, i.generator_point(), i.r);
            RhoReferences {
                matched: rho_reference_for(g, gen, r, seed, repeats, runs, |t, _, s, cap| {
                    rho_reference_negation(g, gen, t, r, s, cap)
                }),
                plain: rho_reference_for(g, gen, r, seed, repeats, runs, |t, _, s, cap| {
                    rho_reference(g, gen, t, r, s, cap)
                }),
                candidates: Vec::new(),
                rule: "negation map (A = 2): the only eligible automorphism of this curve",
            }
        }
        Instance::Binary(i) => {
            let bg = BinaryGroup(&i.fast);
            let (gen, r) = (i.generator, i.r);
            let negation = rho_reference_for(&bg, gen, r, seed, repeats, runs, |t, _, s, cap| {
                rho_reference_negation(&bg, gen, t, r, s, cap)
            });
            let plain = rho_reference_for(&bg, gen, r, seed, repeats, runs, |t, _, s, cap| {
                rho_reference(&bg, gen, t, r, s, cap)
            });
            if i.koblitz.is_none() {
                return RhoReferences {
                    matched: negation,
                    plain,
                    candidates: Vec::new(),
                    rule: "negation map (A = 2): the only eligible automorphism of this curve",
                };
            }
            let signed = rho_reference_for(&bg, gen, r, seed, repeats, runs, |t, planted, s, _| {
                signed_frobenius_rho(i, t, planted, s)
                    .expect("a Koblitz instance carries its curve")
            });
            let mut candidates: Vec<RhoReference> = signed.into_iter().chain(negation).collect();
            // The cheaper of the eligible walks, among those that verified.
            let best = candidates
                .iter()
                .enumerate()
                .filter(|(_, c)| c.all_verified)
                .min_by(|a, b| a.1.mean_s.total_cmp(&b.1.mean_s))
                .map(|(k, _)| k);
            let matched = best.map(|k| candidates.remove(k));
            RhoReferences {
                matched,
                plain,
                candidates,
                rule: "Koblitz curve: the signed-Frobenius walk (A = 2n) and the negation walk (A = 2) on the same seeds; the lower mean S prices the column, the other is kept in rho_reference_candidates",
            }
        }
    }
}

/// One prime-field configuration, over `repeats` targets.
fn run_prime(
    inst: &PrimeInstance,
    spec: &PipelineSpec,
    args: &BenchArgs,
    calib: &Calibration,
    rho_s: Option<f64>,
) -> Result<Vec<RunReport>, String> {
    let fb_name = spec.factor_base.as_str();
    let or_name = spec.oracle.as_str();
    let abscissa = PrimeAbscissaBase { instance: inst };
    let orbit = GlvOrbitBase { instance: inst };
    let base: &dyn FactorBaseBuilder<_> = match fb_name {
        "prime-abscissa" => &abscissa,
        "glv-orbit" => &orbit,
        other => {
            return Err(format!(
                "factor base `{other}` is not available on a prime-field curve; try prime-abscissa or glv-orbit"
            ))
        }
    };
    let mut out = Vec::new();
    for rep in 0..args.repeats.max(1) {
        let mut ops = GroupOps::default();
        let g = inst.generator_point();
        let planted = planted_log(spec.seed, rep, inst.r);
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
        let m = spec.oracle_params.u64_or("m", 2)? as u32;
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(m);
        let oracle: &mut dyn DecompositionOracle<_> = match or_name {
            "subtract" => &mut subtract,
            "mitm" => &mut mitm,
            other => {
                return Err(format!(
                    "oracle `{other}` is not available on a prime-field curve; try subtract or mitm"
                ))
            }
        };
        out.push(run_pipeline(
            &ctx, &spec, base, oracle, planted, calib, rho_s,
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
    rho_s: Option<f64>,
    frozen: &mut Vec<Value>,
) -> Result<Vec<RunReport>, String> {
    let fb_name = spec.factor_base.as_str();
    let or_name = spec.oracle.as_str();
    let g = BinaryGroup(&inst.fast);
    let subspace = BinarySubspaceBase { instance: inst };
    let orbit = KoblitzOrbitBase { instance: inst };
    let trace_zero = KoblitzTraceZeroBase { instance: inst };
    let symmetrised_base = KoblitzSymmetrisedBase { instance: inst };
    let base: &dyn FactorBaseBuilder<BinaryGroup> = match fb_name {
        "binary-subspace" => &subspace,
        "koblitz-orbit" => &orbit,
        "koblitz-trace-zero" => &trace_zero,
        "koblitz-symmetrised" => &symmetrised_base,
        other => {
            return Err(format!(
                "factor base `{other}` is not available on a binary curve; try binary-subspace, koblitz-orbit, koblitz-trace-zero or koblitz-symmetrised"
            ))
        }
    };
    let mut out = Vec::new();
    for rep in 0..args.repeats.max(1) {
        let mut ops = GroupOps::default();
        let planted = planted_log(spec.seed, rep, inst.r);
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
        // The base does not depend on the target, so the first repeat's
        // is the one every repeat ran over.
        if rep == 0 && args.factor_base_out.is_some() {
            let mut doc = frozen_factor_base(inst, &spec, base, &ctx)?;
            doc["configuration"] = json!(spec.label());
            frozen.push(doc);
        }
        let m = spec.oracle_params.u64_or("m", 2)? as u32;
        // The symmetrised oracle rebuilds F_u from the base's divisor and
        // checks it point for point; spare the user typing it twice.
        if or_name == "symmetrised" && spec.oracle_params.get("divisor").is_none() {
            if let Some(d) = spec.factor_base_params.get("divisor") {
                let d = d.to_string();
                spec.oracle_params.set("divisor", d);
            }
        }
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(m);
        let mut frob = FrobeniusMitmOracle::new(m, inst);
        let mut symmetrised = match or_name {
            "symmetrised" => Some(SymmetrisedOracle::new(m, inst)),
            _ => None,
        };
        let mut algebraic = match or_name {
            "descent-algebraic" => {
                let raw = spec.solver.as_deref().ok_or(
                    "descent-algebraic needs --solver (try buchberger-f2, sat-cdcl or exhaustive)",
                )?;
                let (sname, sparams) = parse_plugin(raw)?;
                let budget = (spec.solver_budget_seconds > 0)
                    .then(|| std::time::Duration::from_secs(spec.solver_budget_seconds));
                Some(DescentAlgebraicOracle::new(
                    m,
                    inst,
                    solver_by_name(&sname)?,
                    sparams,
                    budget,
                ))
            }
            _ => None,
        };
        let oracle: &mut dyn DecompositionOracle<BinaryGroup> = match or_name {
            "subtract" => &mut subtract,
            "mitm" => &mut mitm,
            "mitm-frobenius" => &mut frob,
            "descent-algebraic" => algebraic.as_mut().expect("built above"),
            "symmetrised" => symmetrised.as_mut().expect("built above"),
            other => return Err(format!("oracle `{other}` is not available here")),
        };
        out.push(run_pipeline(
            &ctx, &spec, base, oracle, planted, calib, rho_s,
        )?);
    }
    Ok(out)
}

/// The instance a sweep runs every configuration on, built once.
#[allow(clippy::large_enum_variant)]
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
    /// `prime` only: the curve family (`generic`, `j0`, `j1728`, `d7`,
    /// `d8`); the `--family` default when absent.
    #[serde(default)]
    family: Option<String>,
    /// `char2` only: the largest cofactor the random curve may have;
    /// the `--max-cofactor` default when absent.
    #[serde(default)]
    max_cofactor: Option<u64>,
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
) -> Result<PipelineSpec, String> {
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
    Ok(PipelineSpec {
        factor_base: fb_name,
        factor_base_params: fb_params,
        oracle: or_name,
        oracle_params: or_params,
        solver: cfg
            .get("solver")
            .cloned()
            .or_else(|| defaults.solver.clone()),
        solver_params: Params::default(),
        targets: Targets::parse(&tg)?,
        linalg: cfg
            .get("linalg")
            .cloned()
            .unwrap_or_else(|| defaults.linalg.clone()),
        max_trials: defaults.max_trials,
        seed: defaults.seed,
        solver_budget_seconds: defaults.solver_budget_seconds,
    })
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
        let (n, params) = parse_plugin(name)?;
        solver_by_name(&n)?;
        validate_solver_params(&n, &params)?;
    }

    // Either a sweep file or the flags describe the work.
    let (instance_spec, configs, args) = match &args.sweep {
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
            if let Some(v) = file.instance.max_cofactor {
                a.max_cofactor = v;
            }
            if let Some(v) = &file.instance.family {
                a.family = v.clone();
            }
            let mut configs = file.configurations.clone();
            configs.extend(expand_matrix(&file.matrix));
            if configs.is_empty() {
                return Err("the sweep file lists no configurations and no matrix".into());
            }
            (
                (file.instance.regime.clone(), file.instance.degree),
                configs,
                a,
            )
        }
        None => {
            let (regime, degree) = match (args.bits, args.char2_degree, args.koblitz_degree) {
                (Some(b), None, None) => ("prime".to_string(), b),
                (None, Some(n), None) => ("char2".to_string(), n),
                (None, None, Some(n)) => ("koblitz".to_string(), n),
                _ => return Err(
                    "choose exactly one of --bits, --char2-degree, --koblitz-degree, or --sweep"
                        .into(),
                ),
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

    // A sweep may override the CLI default solver per configuration.
    // Reject malformed SAT encoding before calibrations, rho or a target.
    for cfg in &configs {
        if let Some(raw) = cfg.get("solver").or(args.solver.as_ref()) {
            let (name, params) = parse_plugin(raw)?;
            solver_by_name(&name)?;
            validate_solver_params(&name, &params)?;
        }
    }

    let (regime, degree) = instance_spec;
    let instance = match regime.as_str() {
        "prime" => Instance::Prime(match CmFamily::parse(&args.family)? {
            CmFamily::Generic => roster_prime_instance(degree)
                .ok_or_else(|| format!("no prime instance at {degree} bits"))?,
            family => generate_cm_instance(family, degree, args.seed, args.max_cofactor)?,
        }),
        "char2" => Instance::Binary(
            random_binary_instance(degree, args.seed, args.max_cofactor).ok_or_else(|| {
                format!(
                    "no random binary instance at degree {degree} with cofactor at most {}",
                    args.max_cofactor
                )
            })?,
        ),
        "koblitz" => Instance::Binary(match args.koblitz_a {
            Some(a) => koblitz_instance(a, degree)
                .ok_or_else(|| format!("no Koblitz instance K_{a} at degree {degree}"))?,
            None => koblitz_instance(1, degree)
                .or_else(|| koblitz_instance(0, degree))
                .ok_or_else(|| format!("no Koblitz instance at degree {degree}"))?,
        }),
        other => {
            return Err(format!(
                "unknown regime `{other}`; try prime, char2 or koblitz"
            ))
        }
    };

    let (calib, pins) = match &instance {
        Instance::Prime(i) => {
            calibration_for(&i.curve, &sample_prime_points(i), "prime", &i.name, i.r)
        }
        Instance::Binary(i) => calibration_for_binary(i, &regime),
    };

    let started = std::time::Instant::now();
    // The method boundary on this instance, before any configuration:
    // every row's `vs rho` is its S over the matched reference's mean.
    let references = matched_rho(&instance, args.seed, args.repeats, args.rho_runs);
    let rho = &references.matched;
    let rho_s = rho.as_ref().filter(|r| r.all_verified).map(|r| r.mean_s);
    if !json_only {
        if let Some(r) = rho {
            eprintln!(
                "  rho reference (A = {}): S = {:.3} (mean of {}, {:.3}–{:.3}){}",
                r.automorphisms,
                r.mean_s,
                r.runs,
                r.min_s,
                r.max_s,
                if r.all_verified {
                    ""
                } else {
                    " — a run failed to verify; the column is left empty"
                }
            );
        }
        if let Some(p) = &references.plain {
            eprintln!("  before mark, plain walk (A = 1): S = {:.3}", p.mean_s);
        }
        for c in &references.candidates {
            eprintln!("  also run (A = {}): S = {:.3}", c.automorphisms, c.mean_s);
        }
    }
    if args.factor_base_out.is_some() && matches!(instance, Instance::Prime(_)) {
        return Err("--factor-base-out freezes binary and Koblitz bases only".into());
    }
    let mut rows: Vec<RunReport> = Vec::new();
    let mut failures: Vec<Value> = Vec::new();
    let mut frozen: Vec<Value> = Vec::new();
    for cfg in &configs {
        let spec = spec_from(cfg, &args)?;
        if !json_only {
            eprintln!("  {} …", spec.label());
        }
        // A configuration that does not apply to this regime is
        // reported, not fatal: a sweep over a matrix will contain
        // combinations that do not exist, and the useful output is the
        // ones that do plus a note on the ones that do not.
        let outcome = match &instance {
            Instance::Prime(i) => run_prime(i, &spec, &args, &calib, rho_s),
            Instance::Binary(i) => run_binary(i, &spec, &args, &calib, rho_s, &mut frozen),
        };
        match outcome {
            Ok(mut got) => rows.append(&mut got),
            Err(why) => failures.push(json!({"configuration": spec.label(), "why": why})),
        }
    }

    if let Some(path) = &args.factor_base_out {
        let doc = json!({
            "schema_version": 1,
            "operation": "bench-factor-bases",
            "instance": instance.name(),
            "seed": args.seed,
            "factor_bases": frozen,
        });
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
        serde_json::to_writer_pretty(file, &doc).map_err(|e| e.to_string())?;
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
        // The conversion every non-addition unit was priced at, so a
        // frozen report re-derives its own GAE: nanoseconds per unit on
        // this host, with the units the repository's table pinned over
        // the measurement named apart from the ones left measured.
        "calibration": calib,
        "calibration_pins": pins,
        "rho_reference": rho,
        "rho_reference_rule": references.rule,
        "rho_reference_plain": references.plain,
        "rho_reference_candidates": references.candidates,
        "configurations_run": rows.len(),
        "configurations_skipped": failures,
        "all_verified": all_verified,
        "elapsed_seconds": started.elapsed().as_secs_f64(),
        "markdown": markdown,
        "rows": rows,
    }))
}
