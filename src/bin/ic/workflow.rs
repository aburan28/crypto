//! Resumable, parameter-file-driven pipeline: select → collect → logs → solve.
//!
//! The number-field-sieve tools run as a sequence of stages with their
//! outputs on disk, so a run can be stopped, inspected and resumed
//! without redoing finished work.  `ic workflow` gives the Koblitz
//! index calculus the same shape:
//!
//! 1. **select** — choose the factor base: an explicit recipe from the
//!    parameter file, or the best-by-census candidate of the
//!    factor-base search; materialised and written as
//!    `factor_base.json`.
//! 2. **collect** — draw and decompose probes in **work units**, each a
//!    slice of the seed's probe sequence, written as
//!    `relations/unit-NNNNN.json`.  Units are independent: any machine
//!    with the parameter file (and `factor_base.json` for a searched
//!    base) can run `--collect-units 3-7` and drop its files into the
//!    run directory, the way sieving clients feed a number-field-sieve
//!    server.  Every relation carries only the probe scalar and the
//!    factor-base point indices; the driver re-verifies each one in
//!    the group before it is used, so a bad worker cannot poison the
//!    database.
//! 3. **logs** — merge the units, verify and deduplicate the relations,
//!    and precompute the factor-base logarithm database (`logs.json`),
//!    every column certified by `[x]G == R`.  When the relations do not
//!    yet determine every column the driver collects further units up
//!    to `collection.max_units`.
//! 4. **solve** — descend each target with one relation reusing the
//!    database, appending to `solutions.json` after every target so an
//!    interrupted run resumes at the first unsolved one.
//!
//! `state.json` in the run directory records the parameter digest and
//! each stage's status.  A rerun in the same directory loads whatever
//! artifacts exist, re-verifies them against the reconstructed curve
//! (a stale or tampered artifact is rejected, never trusted), and
//! continues from the first stage that is not complete.  A parameter
//! file that no longer matches the recorded digest is refused, so one
//! directory never silently mixes two experiments.  Artifacts are
//! written atomically (temp file, then rename).
//!
//! Every target is a synthetic known-answer instance built from the
//! parameter file; imported points are never solved.

use super::experiment::{
    self, log_table_from_doc, log_table_to_doc, FactorBaseDocument, LinearAlgebraMode,
    LogTableDocument, Solver,
};
use crypto_lib::cryptanalysis::koblitz_sparse_la::SparseSolveOptions;
use clap::{Args, ValueEnum};
use crypto_lib::cryptanalysis::koblitz_factor_base_search::{
    search, Candidate, FactorBaseSpec, Family, SearchOptions,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    individual_log_with_pair_table, solve_factor_base_logs_from_relations, CollectedRelation,
    DecompositionStrategy, FrobeniusFactorBase, KoblitzCurve, KoblitzIcOptions, PairSumTable,
    RelationCollector, RelationWorkUnit,
};
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{Read, Write},
    path::{Path, PathBuf},
    time::Instant,
};

// ── Parameter file ─────────────────────────────────────────────────

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CurveParams {
    pub degree: u32,
    #[serde(default)]
    pub curve_a: u8,
    /// Subfield degree `k` (`q = 2^k`); 1 is the Koblitz family.
    #[serde(default = "one_u32", skip_serializing_if = "is_one_u32")]
    pub subfield: u32,
    /// Coordinate of `b` in the subfield basis (1 for Koblitz curves).
    #[serde(default = "one_u64", skip_serializing_if = "is_one_u64")]
    pub curve_b: u64,
}
fn one_u32() -> u32 {
    1
}
fn one_u64() -> u64 {
    1
}
fn is_one_u32(v: &u32) -> bool {
    *v == 1
}
fn is_one_u64(v: &u64) -> bool {
    *v == 1
}
/// Whether a document's curve fields name `c`.
fn same_curve(degree: u32, curve_a: u8, subfield: u32, curve_b: u64, c: &KoblitzCurve) -> bool {
    degree == c.n && curve_a == c.a && subfield == c.k && curve_b == c.b_index
}

/// How the select stage obtains the factor base.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub enum FactorBaseSource {
    /// Use this recipe as is.
    Spec { spec: FactorBaseSpec },
    /// Run the factor-base search and take its best-by-census candidate.
    Search {
        #[serde(default = "default_family")]
        family: String,
        #[serde(default = "default_min_dimension")]
        min_dimension: u32,
        #[serde(default = "default_max_dimension")]
        max_dimension: u32,
        #[serde(default = "default_max_abscissae")]
        max_abscissae: usize,
        #[serde(default = "default_union_min_seed")]
        union_min_seed: u32,
        #[serde(default = "default_union_max_seed")]
        union_max_seed: u32,
        #[serde(default = "default_union_samples")]
        union_samples: usize,
        #[serde(default = "default_targets")]
        targets: usize,
        #[serde(default = "default_exhaustive_cap")]
        exhaustive_cap: u64,
        #[serde(default = "default_extra_relations")]
        extra_relations: usize,
        #[serde(default = "default_true")]
        prune: bool,
        #[serde(default)]
        saturate: bool,
    },
}
fn default_family() -> String {
    "all".into()
}
fn default_min_dimension() -> u32 {
    3
}
fn default_max_dimension() -> u32 {
    10
}
fn default_max_abscissae() -> usize {
    2048
}
fn default_union_min_seed() -> u32 {
    2
}
fn default_union_max_seed() -> u32 {
    5
}
fn default_union_samples() -> usize {
    3
}
fn default_targets() -> usize {
    1024
}
fn default_exhaustive_cap() -> u64 {
    4096
}
fn default_extra_relations() -> usize {
    2
}
fn default_true() -> bool {
    true
}

/// One synthetic known-answer target: an explicit scalar, or one drawn
/// reproducibly from a seed.  Exactly one field must be set.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TargetSpec {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub known_log: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub random_seed: Option<u64>,
}

/// The linear algebra of the logs stage: `dense` elimination or
/// `sparse` filtering + block Wiedemann with its knobs.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LinearAlgebraParams {
    pub mode: LinearAlgebraMode,
    pub sparse: SparseSolveOptions,
}

/// Relation collection in work units.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct CollectionParams {
    /// Probes per work unit.
    pub unit_trials: u64,
    /// Units the driver collects before the first solve attempt.
    pub units: usize,
    /// Units the driver may collect in total when the relations do not
    /// yet determine every column.
    pub max_units: usize,
}
impl Default for CollectionParams {
    fn default() -> Self {
        Self {
            unit_trials: 4096,
            units: 4,
            max_units: 64,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkflowParams {
    pub schema_version: u32,
    #[serde(default)]
    pub name: String,
    pub curve: CurveParams,
    #[serde(default = "default_summands")]
    pub summands: u8,
    #[serde(default = "default_solver")]
    pub solver: Solver,
    #[serde(default = "default_seed")]
    pub seed: u64,
    #[serde(default = "default_max_trials")]
    pub max_trials: u32,
    #[serde(default)]
    pub linear_algebra: LinearAlgebraParams,
    #[serde(default)]
    pub collection: CollectionParams,
    pub factor_base: FactorBaseSource,
    #[serde(default)]
    pub targets: Vec<TargetSpec>,
}
fn default_summands() -> u8 {
    2
}
fn default_solver() -> Solver {
    Solver::PairTable
}
fn default_seed() -> u64 {
    0x4b_6f_62_6c_69_74_7a_00
}
fn default_max_trials() -> u32 {
    200_000
}

pub fn load_params(path: &Path) -> Result<WorkflowParams, String> {
    let mut data = Vec::new();
    File::open(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .take(1_048_577)
        .read_to_end(&mut data)
        .map_err(|e| e.to_string())?;
    if data.len() > 1_048_576 {
        return Err("workflow parameter file exceeds 1 MiB".into());
    }
    let p: WorkflowParams =
        serde_json::from_slice(&data).map_err(|e| format!("invalid workflow parameters: {e}"))?;
    if p.schema_version != 1 {
        return Err("unsupported workflow schema version".into());
    }
    if p.summands < 2 || p.summands > 4 {
        return Err("summands must be 2, 3 or 4".into());
    }
    if p.max_trials == 0 || p.max_trials > 1_000_000 {
        return Err("max_trials must be 1..=1000000".into());
    }
    experiment::field_degree(&p.curve.degree.to_string())?;
    if p.curve.subfield == 1 && p.curve.curve_a > 1 {
        return Err("curve_a must be 0 or 1 for a Koblitz curve".into());
    }
    let col = &p.collection;
    if col.unit_trials == 0 || col.unit_trials > 100_000_000 {
        return Err("collection.unit_trials must be 1..=100000000".into());
    }
    if col.units == 0 || col.max_units < col.units || col.max_units > 100_000 {
        return Err("collection.units must be at least 1 and at most collection.max_units (≤ 100000)".into());
    }
    for (i, t) in p.targets.iter().enumerate() {
        if t.known_log.is_some() == t.random_seed.is_some() {
            return Err(format!(
                "target {i}: set exactly one of known_log or random_seed"
            ));
        }
    }
    Ok(p)
}

/// Digest of the canonical parameter JSON; one run directory serves one
/// parameter set.
fn params_digest(p: &WorkflowParams) -> String {
    let canonical = serde_json::to_vec(p).expect("parameters serialise");
    blake3::hash(&canonical).to_hex().to_string()
}

// ── Persistent state ───────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum Stage {
    Select,
    Collect,
    Logs,
    Solve,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum StageStatus {
    #[default]
    Pending,
    Complete,
    Failed,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StageState {
    #[serde(default)]
    pub status: StageStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub artifact: Option<String>,
    #[serde(default)]
    pub elapsed_seconds: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkflowState {
    pub schema_version: u32,
    pub params_digest: String,
    pub name: String,
    #[serde(default)]
    pub select: StageState,
    #[serde(default)]
    pub collect: StageState,
    #[serde(default)]
    pub logs: StageState,
    #[serde(default)]
    pub solve: StageState,
    #[serde(default)]
    pub solved_targets: usize,
    #[serde(default)]
    pub total_targets: usize,
    #[serde(default)]
    pub runs: u32,
    /// Valid work-unit files present at the last run.
    #[serde(default)]
    pub units_collected: usize,
    /// Relations in those units (before verification and deduplication).
    #[serde(default)]
    pub relations_collected: usize,
}

/// One work unit's relations, bound to the run and the factor base so a
/// file from another experiment is ignored rather than merged.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelationUnitDocument {
    pub schema_version: u32,
    pub params_digest: String,
    pub degree: u32,
    pub curve_a: u8,
    #[serde(default = "one_u32", skip_serializing_if = "is_one_u32")]
    pub subfield: u32,
    #[serde(default = "one_u64", skip_serializing_if = "is_one_u64")]
    pub curve_b: u64,
    pub spec: FactorBaseSpec,
    pub summands: u8,
    pub unit: usize,
    pub seed: u64,
    pub start: u64,
    pub count: u64,
    pub relations: Vec<CollectedRelation>,
    pub elapsed_seconds: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Solution {
    pub index: usize,
    pub expected: String,
    pub recovered: Option<String>,
    pub verified: bool,
    pub descent_trials: usize,
    pub elapsed_seconds: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SolutionsDocument {
    pub schema_version: u32,
    pub degree: u32,
    pub curve_a: u8,
    #[serde(default = "one_u32", skip_serializing_if = "is_one_u32")]
    pub subfield: u32,
    #[serde(default = "one_u64", skip_serializing_if = "is_one_u64")]
    pub curve_b: u64,
    pub params_digest: String,
    pub solutions: Vec<Solution>,
}

const STATE_FILE: &str = "state.json";
const FACTOR_BASE_FILE: &str = "factor_base.json";
const LOGS_FILE: &str = "logs.json";
const SOLUTIONS_FILE: &str = "solutions.json";
const RELATIONS_DIR: &str = "relations";

fn unit_path(dir: &Path, unit: usize) -> PathBuf {
    dir.join(RELATIONS_DIR).join(format!("unit-{unit:05}.json"))
}

/// A worker's list of work units, as given on the command line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnitList(pub Vec<usize>);

/// Parse a worker's unit list: `3`, `0-3`, `0-3,7,9-10`.
pub fn parse_units(text: &str) -> Result<UnitList, String> {
    let mut out = Vec::new();
    for part in text.split(',').map(str::trim).filter(|p| !p.is_empty()) {
        let (lo, hi) = match part.split_once('-') {
            Some((a, b)) => (a.trim().parse::<usize>(), b.trim().parse::<usize>()),
            None => (part.parse::<usize>(), part.parse::<usize>()),
        };
        let (lo, hi) = (lo.map_err(|e| format!("{part:?}: {e}"))?, hi.map_err(|e| format!("{part:?}: {e}"))?);
        if lo > hi {
            return Err(format!("{part:?}: empty range"));
        }
        if hi - lo > 100_000 {
            return Err(format!("{part:?}: range too wide"));
        }
        out.extend(lo..=hi);
    }
    out.sort_unstable();
    out.dedup();
    if out.is_empty() {
        return Err("no work units given".into());
    }
    Ok(UnitList(out))
}

/// Write a JSON document atomically: temp file in the same directory,
/// then rename over the destination.
fn write_atomic<T: Serialize>(path: &Path, value: &T) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|e| e.to_string())?;
    let tmp = path.with_extension("json.tmp");
    {
        let mut f = File::create(&tmp).map_err(|e| format!("{}: {e}", tmp.display()))?;
        f.write_all(text.as_bytes()).map_err(|e| e.to_string())?;
        f.write_all(b"\n").map_err(|e| e.to_string())?;
        f.sync_all().map_err(|e| e.to_string())?;
    }
    fs::rename(&tmp, path).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_json<T: for<'de> Deserialize<'de>>(path: &Path) -> Result<T, String> {
    let mut data = Vec::new();
    File::open(path)
        .map_err(|e| format!("{}: {e}", path.display()))?
        .take(16_777_217)
        .read_to_end(&mut data)
        .map_err(|e| e.to_string())?;
    if data.len() > 16_777_216 {
        return Err(format!("{} exceeds 16 MiB", path.display()));
    }
    serde_json::from_slice(&data).map_err(|e| format!("{}: {e}", path.display()))
}

// ── Driver ─────────────────────────────────────────────────────────

#[derive(Clone, Debug, Args)]
pub struct WorkflowArgs {
    /// Workflow parameter file (JSON, schema_version 1).
    #[arg(long)]
    pub params: PathBuf,
    /// Run directory holding state and artifacts; created if missing.
    #[arg(long)]
    pub dir: PathBuf,
    /// Stop after this stage (useful for staged or interrupted runs).
    #[arg(long, value_enum)]
    pub stop_after: Option<Stage>,
    /// Act as a collection worker: run only these work units (`3`,
    /// `0-3`, `0-3,7`), write their relation files, and stop.
    #[arg(long, value_parser = parse_units)]
    pub collect_units: Option<UnitList>,
}

fn known_scalar(c: &KoblitzCurve, t: &TargetSpec) -> Result<BigUint, String> {
    let k = if let Some(s) = &t.known_log {
        super::params::number(s)?
    } else {
        let seed = t.random_seed.expect("validated: one of known_log/random_seed");
        let mut rng = StdRng::seed_from_u64(seed ^ 0x534f_4c56_4552_5447);
        BigUint::from(rng.gen_range(1..c.subgroup_order.to_u64_digits()[0]))
    };
    if k.is_zero_() || k >= c.subgroup_order {
        return Err(format!(
            "target scalar must be nonzero and smaller than the subgroup order {}",
            c.subgroup_order
        ));
    }
    Ok(k)
}

trait IsZeroExt {
    fn is_zero_(&self) -> bool;
}
impl IsZeroExt for BigUint {
    fn is_zero_(&self) -> bool {
        use num_traits::Zero;
        self.is_zero()
    }
}

fn search_options(src: &FactorBaseSource, p: &WorkflowParams) -> Result<SearchOptions, String> {
    let FactorBaseSource::Search {
        family,
        min_dimension,
        max_dimension,
        max_abscissae,
        union_min_seed,
        union_max_seed,
        union_samples,
        targets,
        exhaustive_cap,
        extra_relations,
        prune,
        saturate,
    } = src
    else {
        return Err("not a search source".into());
    };
    let families = match family.as_str() {
        "factor" => vec![Family::Factor],
        "divisor" => vec![Family::Divisor],
        "union" => vec![Family::Union],
        "all" => vec![Family::Factor, Family::Divisor, Family::Union],
        other => return Err(format!("unknown factor-base family {other:?}")),
    };
    if min_dimension > max_dimension || union_min_seed > union_max_seed {
        return Err("dimension windows must be non-decreasing".into());
    }
    if *max_abscissae > experiment::MAX_ABSCISSAE {
        return Err(format!(
            "max_abscissae exceeds the materialization limit {}",
            experiment::MAX_ABSCISSAE
        ));
    }
    Ok(SearchOptions {
        m: p.summands as usize,
        min_dimension: *min_dimension,
        max_dimension: *max_dimension,
        max_abscissae: *max_abscissae,
        families,
        union_seed_dimensions: (*union_min_seed, *union_max_seed),
        union_samples: *union_samples,
        sample_targets: *targets,
        exhaustive_cap: *exhaustive_cap,
        extra_relations: *extra_relations,
        prune: *prune,
        saturate: *saturate,
        projected_columns: true,
        seed: p.seed,
    })
}

fn candidate_json(c: &Candidate) -> Value {
    let census = c.census.as_ref();
    json!({"spec":c.spec,"family":c.family,"points":c.points,"columns":c.unknowns,
        "coverage":census.map(|x| x.coverage),
        "expected_trials":census.map(|x| if x.expected_trials.is_finite(){json!(x.expected_trials)}else{Value::Null})})
}

/// Run (or resume) the workflow described by `args.params` in `args.dir`.
pub fn run(args: WorkflowArgs, quiet: bool) -> Result<Value, String> {
    let begin = Instant::now();
    let p = load_params(&args.params)?;
    let digest = params_digest(&p);
    fs::create_dir_all(&args.dir).map_err(|e| format!("{}: {e}", args.dir.display()))?;
    let state_path = args.dir.join(STATE_FILE);
    let fb_path = args.dir.join(FACTOR_BASE_FILE);
    let logs_path = args.dir.join(LOGS_FILE);
    let sol_path = args.dir.join(SOLUTIONS_FILE);

    // Load or initialise state; refuse to mix parameter sets.
    let mut state: WorkflowState = if state_path.exists() {
        let s: WorkflowState = read_json(&state_path)?;
        if s.schema_version != 1 {
            return Err("unsupported workflow state schema".into());
        }
        if s.params_digest != digest {
            return Err(format!(
                "run directory {} belongs to a different parameter set (digest {}); use a fresh directory",
                args.dir.display(),
                &s.params_digest[..16]
            ));
        }
        s
    } else {
        WorkflowState {
            schema_version: 1,
            params_digest: digest.clone(),
            name: p.name.clone(),
            select: StageState::default(),
            collect: StageState::default(),
            logs: StageState::default(),
            solve: StageState::default(),
            solved_targets: 0,
            total_targets: p.targets.len(),
            runs: 0,
            units_collected: 0,
            relations_collected: 0,
        }
    };
    state.runs += 1;
    state.total_targets = p.targets.len();
    let resumed = state.runs > 1;
    let say = |msg: &str| {
        if !quiet {
            println!("{msg}");
            let _ = std::io::stdout().flush();
        }
    };
    say(&format!(
        "ic — workflow {:?} in {} ({}); {}; {} summands; engine {}",
        p.name,
        args.dir.display(),
        if resumed { format!("resuming, run {}", state.runs) } else { "fresh".into() },
        experiment::curve_label(p.curve.degree, p.curve.curve_a, p.curve.subfield, p.curve.curve_b),
        p.summands,
        p.solver.name()
    ));

    let c = experiment::curve(p.curve.degree, p.curve.curve_a, p.curve.subfield, p.curve.curve_b)?;
    let mut stage_reports: Vec<Value> = Vec::new();
    let mut overall_failed: Option<String> = None;

    // ── Stage 1: select ────────────────────────────────────────────
    let t0 = Instant::now();
    let mut select_ran = false;
    let (spec, fb): (FactorBaseSpec, FrobeniusFactorBase) =
        if state.select.status == StageStatus::Complete && fb_path.exists() {
            let doc: FactorBaseDocument = read_json(&fb_path)?;
            if !doc.matches(&c) {
                return Err("factor_base.json does not belong to this curve".into());
            }
            let fb = experiment::materialize(&c, &doc.spec)?;
            say(&format!("[1/4] select: reused {}", FACTOR_BASE_FILE));
            (doc.spec, fb)
        } else {
            select_ran = true;
            let (spec, search_report) = match &p.factor_base {
                FactorBaseSource::Spec { spec } => (spec.clone(), None),
                src @ FactorBaseSource::Search { .. } => {
                    let opts = search_options(src, &p)?;
                    say("[1/4] select: running factor-base search …");
                    let report = search(&c, &opts);
                    let best = report.best().ok_or(
                        "factor-base search found no candidate that decomposes any target",
                    )?;
                    let top: Vec<Value> = report.candidates.iter().take(5).map(candidate_json).collect();
                    (best.spec.clone(), Some(json!({"candidates_scored":report.candidates.len(),
                        "targets":report.targets,"exhaustive_targets":report.exhaustive_targets,
                        "elapsed_ms":report.elapsed_ms,"top":top})))
                }
            };
            let fb = experiment::materialize(&c, &spec)?;
            let doc = FactorBaseDocument {
                schema_version: 1,
                degree: c.n,
                curve_a: c.a,
                subfield: c.k,
                curve_b: c.b_index,
                spec: spec.clone(),
            };
            write_atomic(&fb_path, &doc)?;
            state.select = StageState {
                status: StageStatus::Complete,
                artifact: Some(FACTOR_BASE_FILE.into()),
                elapsed_seconds: t0.elapsed().as_secs_f64(),
                reason: None,
            };
            write_atomic(&state_path, &state)?;
            if let Some(sr) = search_report {
                stage_reports.push(json!({"stage":"select","status":"complete","ran":true,"search":sr}));
            }
            say(&format!(
                "[1/4] select: complete — {} ({} points, {} signed orbits)",
                serde_json::to_string(&spec).unwrap_or_default(),
                fb.points.len(),
                fb.unknowns()
            ));
            (spec, fb)
        };
    if !select_ran {
        stage_reports.push(json!({"stage":"select","status":"complete","ran":false}));
    } else if stage_reports.last().map_or(true, |v| v["stage"] != "select") {
        stage_reports.push(json!({"stage":"select","status":"complete","ran":true}));
    }
    let columns = crypto_lib::cryptanalysis::koblitz_index_calculus::projected_signed_orbit_count(&c, &fb);
    let factor_base_summary = experiment::factor_base_json(&spec, &fb, columns);
    if args.stop_after == Some(Stage::Select) {
        return Ok(finish(&p, &state, &args, stage_reports, factor_base_summary, None, None, begin, "stopped"));
    }

    // ── Stage 2: collect (work units) ──────────────────────────────
    let t1 = Instant::now();
    let ic = experiment::with_linear_algebra(
        experiment::ic_options(p.solver, p.summands, p.max_trials, p.seed),
        p.linear_algebra.mode,
        p.linear_algebra.sparse,
    );
    let rel_dir = args.dir.join(RELATIONS_DIR);
    fs::create_dir_all(&rel_dir).map_err(|e| format!("{}: {e}", rel_dir.display()))?;
    // The pair table is built at most once per process and shared by
    // collection and descent.
    let mut pair: Option<PairSumTable> = None;
    let (mut units, ignored) = load_units(&rel_dir, &digest, &c, &spec, p.summands)?;
    let worker_mode = args.collect_units.is_some();
    let logs_done = state.logs.status == StageStatus::Complete && logs_path.exists();
    let wanted: Vec<usize> = match args.collect_units.as_ref().map(|u| &u.0) {
        Some(list) => {
            if let Some(&u) = list.iter().find(|&&u| u >= p.collection.max_units) {
                return Err(format!(
                    "work unit {u} is beyond collection.max_units = {}",
                    p.collection.max_units
                ));
            }
            list.iter().copied().filter(|u| !units.contains_key(u)).collect()
        }
        None if logs_done => Vec::new(),
        None => (0..p.collection.units).filter(|u| !units.contains_key(u)).collect(),
    };
    if ignored > 0 {
        say(&format!(
            "[2/4] collect: ignored {ignored} relation file(s) from another run or factor base"
        ));
    }
    let mut units_ran = 0usize;
    let mut relations_now = 0usize;
    let mut trials_now = 0u64;
    if !wanted.is_empty() {
        say(&format!(
            "[2/4] collect: running {} work unit(s) of {} probes ({} already present) …",
            wanted.len(),
            p.collection.unit_trials,
            units.len()
        ));
    }
    for &u in &wanted {
        let doc = collect_unit(&c, &fb, &ic, &mut pair, &p, &digest, &spec, &rel_dir, u)?;
        say(&format!(
            "      unit {u:05}: {} relations from {} probes ({:.2}s)",
            doc.relations.len(),
            doc.count,
            doc.elapsed_seconds
        ));
        units_ran += 1;
        relations_now += doc.relations.len();
        trials_now += doc.count;
        units.insert(u, doc);
        state.units_collected = units.len();
        state.relations_collected = units.values().map(|d| d.relations.len()).sum();
        write_atomic(&state_path, &state)?;
    }
    state.units_collected = units.len();
    state.relations_collected = units.values().map(|d| d.relations.len()).sum();
    let requested_present = (0..p.collection.units).all(|u| units.contains_key(&u));
    state.collect = StageState {
        status: if requested_present || logs_done { StageStatus::Complete } else { StageStatus::Pending },
        artifact: Some(RELATIONS_DIR.into()),
        elapsed_seconds: state.collect.elapsed_seconds + t1.elapsed().as_secs_f64(),
        reason: None,
    };
    write_atomic(&state_path, &state)?;
    let trials_total: u64 = units.values().map(|d| d.count).sum();
    stage_reports.push(json!({"stage":"collect",
        "status":if state.collect.status == StageStatus::Complete {"complete"} else {"partial"},
        "ran":units_ran > 0,"worker":worker_mode,
        "units":{"present":units.len(),"ran_now":units_ran,"requested":p.collection.units,
                 "max":p.collection.max_units,"ignored":ignored,"trials_per_unit":p.collection.unit_trials,
                 "indices":units.keys().copied().collect::<Vec<_>>()},
        "trials_now":trials_now,"relations_now":relations_now,
        "trials_total":trials_total,"relations_total":state.relations_collected,
        "elapsed_seconds":t1.elapsed().as_secs_f64()}));
    if units_ran == 0 {
        say(&format!(
            "[2/4] collect: reused {} unit(s), {} relations from {} probes",
            units.len(),
            state.relations_collected,
            trials_total
        ));
    } else {
        say(&format!(
            "[2/4] collect: {} unit(s) present, {} relations from {} probes ({:.1}s this run)",
            units.len(),
            state.relations_collected,
            trials_total,
            t1.elapsed().as_secs_f64()
        ));
    }
    if worker_mode || args.stop_after == Some(Stage::Collect) {
        return Ok(finish(&p, &state, &args, stage_reports, factor_base_summary, None, None, begin, "stopped"));
    }

    // ── Stage 3: logs (merge, verify, solve) ───────────────────────
    let t2 = Instant::now();
    let table = if logs_done {
        let doc: LogTableDocument = read_json(&logs_path)?;
        let table = log_table_from_doc(&c, &doc)?;
        if doc.spec != spec {
            return Err("logs.json was built over a different factor base".into());
        }
        say(&format!("[3/4] logs: reused {} ({} columns, re-verified)", LOGS_FILE, table.len()));
        stage_reports.push(json!({"stage":"logs","status":"complete","ran":false,"columns":table.len()}));
        Some(table)
    } else {
        say(&format!(
            "[3/4] logs: verifying {} relations from {} unit(s) and solving …",
            state.relations_collected,
            units.len()
        ));
        let mut merged: Vec<CollectedRelation> =
            units.values().flat_map(|d| d.relations.iter().cloned()).collect();
        let mut outcome = solve_factor_base_logs_from_relations(&c, &fb, &ic, &merged);
        let mut extended = 0usize;
        // Undetermined: collect further units up to the budget, solving
        // again after each.
        while let Some((_, report)) = &outcome {
            if report.verified {
                break;
            }
            let next = units.keys().max().map_or(0, |m| m + 1);
            if next >= p.collection.max_units {
                break;
            }
            say(&format!(
                "      {} accepted relations ({} rejected, {} duplicates) do not determine all {} columns; collecting unit {next:05} …",
                report.relations, report.rejected_relations, report.duplicate_relations, report.columns
            ));
            let doc = collect_unit(&c, &fb, &ic, &mut pair, &p, &digest, &spec, &rel_dir, next)?;
            say(&format!(
                "      unit {next:05}: {} relations from {} probes ({:.2}s)",
                doc.relations.len(),
                doc.count,
                doc.elapsed_seconds
            ));
            merged.extend(doc.relations.iter().cloned());
            units.insert(next, doc);
            extended += 1;
            state.units_collected = units.len();
            state.relations_collected = units.values().map(|d| d.relations.len()).sum();
            write_atomic(&state_path, &state)?;
            outcome = solve_factor_base_logs_from_relations(&c, &fb, &ic, &merged);
        }
        let trials_total: u64 = units.values().map(|d| d.count).sum();
        match outcome {
            Some((table, report)) if report.verified => {
                let doc = log_table_to_doc(&c, &spec, p.summands, p.solver, &table);
                write_atomic(&logs_path, &doc)?;
                state.logs = StageState {
                    status: StageStatus::Complete,
                    artifact: Some(LOGS_FILE.into()),
                    elapsed_seconds: t2.elapsed().as_secs_f64(),
                    reason: None,
                };
                write_atomic(&state_path, &state)?;
                stage_reports.push(json!({"stage":"logs","status":"complete","ran":true,
                    "columns":report.columns,"trials":trials_total,"relations":report.relations,
                    "relations_loaded":merged.len(),"rejected":report.rejected_relations,
                    "duplicates":report.duplicate_relations,"units_used":units.len(),"units_extended":extended,
                    "verification_seconds":report.collection_seconds,
                    "linear_algebra":experiment::linear_algebra_json(&report),
                    "elapsed_seconds":t2.elapsed().as_secs_f64()}));
                say(&format!(
                    "[3/4] logs: complete — {} columns certified from {} relations ({} rejected, {} duplicates) out of {} units / {} probes ({:.1}s); {}",
                    report.columns, report.relations, report.rejected_relations, report.duplicate_relations,
                    units.len(), trials_total, t2.elapsed().as_secs_f64(),
                    super::linear_algebra_summary(&experiment::linear_algebra_json(&report))
                ));
                Some(table)
            }
            Some((_, report)) => {
                let reason = format!(
                    "relations from {} units ({} probes; {} accepted, {} rejected, {} duplicates) did not determine every one of {} columns; raise collection.max_units or collection.unit_trials",
                    units.len(), trials_total, report.relations, report.rejected_relations,
                    report.duplicate_relations, report.columns
                );
                state.logs = StageState {
                    status: StageStatus::Failed,
                    artifact: None,
                    elapsed_seconds: t2.elapsed().as_secs_f64(),
                    reason: Some(reason.clone()),
                };
                write_atomic(&state_path, &state)?;
                stage_reports.push(json!({"stage":"logs","status":"failed","ran":true,"reason":reason,
                    "relations_loaded":merged.len(),"rejected":report.rejected_relations,
                    "duplicates":report.duplicate_relations,"units_used":units.len(),"units_extended":extended,
                    "linear_algebra":experiment::linear_algebra_json(&report)}));
                overall_failed = Some(reason);
                None
            }
            None => {
                let reason = "factor base has no usable projected columns for this summand count".to_string();
                state.logs = StageState {
                    status: StageStatus::Failed,
                    artifact: None,
                    elapsed_seconds: t2.elapsed().as_secs_f64(),
                    reason: Some(reason.clone()),
                };
                write_atomic(&state_path, &state)?;
                stage_reports.push(json!({"stage":"logs","status":"failed","ran":true,"reason":reason}));
                overall_failed = Some(reason);
                None
            }
        }
    };
    let Some(table) = table else {
        return Ok(finish(&p, &state, &args, stage_reports, factor_base_summary, None, overall_failed, begin, "failed"));
    };
    if args.stop_after == Some(Stage::Logs) {
        return Ok(finish(&p, &state, &args, stage_reports, factor_base_summary, None, None, begin, "stopped"));
    }

    // ── Stage 3: solve (per-target resumable) ──────────────────────
    let t2 = Instant::now();
    let mut solutions: SolutionsDocument = if sol_path.exists() {
        let d: SolutionsDocument = read_json(&sol_path)?;
        if d.params_digest != digest || !same_curve(d.degree, d.curve_a, d.subfield, d.curve_b, &c) {
            return Err("solutions.json belongs to a different run".into());
        }
        d
    } else {
        SolutionsDocument {
            schema_version: 1,
            degree: c.n,
            curve_a: c.a,
            subfield: c.k,
            curve_b: c.b_index,
            params_digest: digest.clone(),
            solutions: Vec::new(),
        }
    };
    let already: std::collections::HashSet<usize> =
        solutions.solutions.iter().filter(|s| s.verified).map(|s| s.index).collect();
    let pending: Vec<usize> = (0..p.targets.len()).filter(|i| !already.contains(i)).collect();
    say(&format!(
        "[4/4] solve: {} targets, {} already solved, {} pending",
        p.targets.len(),
        already.len(),
        pending.len()
    ));
    // The pair table is shared with collection; built once per process.
    if ic.strategy == DecompositionStrategy::PairTable && !pending.is_empty() && pair.is_none() {
        pair = Some(PairSumTable::build(&c, &fb).ok_or("field too wide for the pair table")?);
    }
    let mut solved_now = 0usize;
    let mut failed_now = 0usize;
    for i in pending {
        let t = Instant::now();
        let k = known_scalar(&c, &p.targets[i])?;
        let q = c.mul(c.generator(), &k);
        let outcome = individual_log_with_pair_table(&c, &fb, &table, &q, &ic, pair.as_ref());
        let (recovered, trials) = match outcome {
            Some((d, r)) => (Some(d), r.trials),
            None => (None, 0),
        };
        let verified = recovered.as_ref() == Some(&k)
            && recovered.as_ref().is_some_and(|d| c.mul(c.generator(), d) == q);
        // Replace any earlier unverified attempt for this index.
        solutions.solutions.retain(|s| s.index != i);
        solutions.solutions.push(Solution {
            index: i,
            expected: k.to_string(),
            recovered: recovered.map(|d| d.to_string()),
            verified,
            descent_trials: trials,
            elapsed_seconds: t.elapsed().as_secs_f64(),
        });
        solutions.solutions.sort_by_key(|s| s.index);
        if verified {
            solved_now += 1;
        } else {
            failed_now += 1;
        }
        state.solved_targets = solutions.solutions.iter().filter(|s| s.verified).count();
        write_atomic(&sol_path, &solutions)?;
        write_atomic(&state_path, &state)?;
        say(&format!(
            "      target {i}: {} ({} descent trials, {:.2}s)",
            if verified { "verified" } else { "FAILED" },
            trials,
            t.elapsed().as_secs_f64()
        ));
    }
    let all_verified = state.solved_targets == p.targets.len();
    state.solve = StageState {
        status: if all_verified { StageStatus::Complete } else { StageStatus::Failed },
        artifact: Some(SOLUTIONS_FILE.into()),
        elapsed_seconds: t2.elapsed().as_secs_f64(),
        reason: (!all_verified).then(|| format!("{} of {} targets unsolved", p.targets.len() - state.solved_targets, p.targets.len())),
    };
    write_atomic(&state_path, &state)?;
    stage_reports.push(json!({"stage":"solve","status":if all_verified{"complete"}else{"failed"},"ran":true,
        "targets":p.targets.len(),"already_solved":already.len(),"solved_now":solved_now,"failed_now":failed_now,
        "elapsed_seconds":t2.elapsed().as_secs_f64()}));
    let status = if all_verified { "complete" } else { "failed" };
    Ok(finish(&p, &state, &args, stage_reports, factor_base_summary, Some(&solutions), overall_failed, begin, status))
}

#[allow(clippy::too_many_arguments)]
fn finish(
    p: &WorkflowParams,
    state: &WorkflowState,
    args: &WorkflowArgs,
    stages: Vec<Value>,
    factor_base: Value,
    solutions: Option<&SolutionsDocument>,
    failure: Option<String>,
    begin: Instant,
    status: &str,
) -> Value {
    json!({"schema_version":1,"operation":"workflow","status":status,
        "evidence_scope":"synthetic_known_answer",
        "name":p.name,"degree":p.curve.degree,"curve_a":p.curve.curve_a,"subfield":p.curve.subfield,"curve_b":p.curve.curve_b,"summands":p.summands,"solver":p.solver,
        "params_digest":state.params_digest,"run_directory":args.dir.display().to_string(),"run_number":state.runs,
        "resumed":state.runs>1,"stop_after":args.stop_after,
        "factor_base":factor_base,
        "stages":stages,
        "state":state,
        "solutions":solutions.map(|s| json!({"count":s.solutions.len(),
            "verified":s.solutions.iter().filter(|x| x.verified).count(),
            "items":s.solutions})),
        "failure":failure,
        "artifacts":{"state":STATE_FILE,"factor_base":FACTOR_BASE_FILE,"relations":RELATIONS_DIR,"logs":LOGS_FILE,"solutions":SOLUTIONS_FILE},
        "elapsed_seconds":begin.elapsed().as_secs_f64(),"resources":experiment::resources(),
        "scope":"resumable staged pipeline: select → collect → logs → solve; completed stages, work units and solved targets are reused on rerun after re-verification against the reconstructed curve",
        "limitations":["No imported target was used.","This run does not establish scaling or challenge readiness."]})
}

/// Read every valid work-unit file of this run; files belonging to
/// another parameter set, curve or factor base — or unreadable ones —
/// are counted and ignored.
fn load_units(
    rel_dir: &Path,
    digest: &str,
    c: &KoblitzCurve,
    spec: &FactorBaseSpec,
    summands: u8,
) -> Result<(BTreeMap<usize, RelationUnitDocument>, usize), String> {
    let mut units = BTreeMap::new();
    let mut ignored = 0usize;
    let entries = match fs::read_dir(rel_dir) {
        Ok(e) => e,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok((units, 0)),
        Err(e) => return Err(format!("{}: {e}", rel_dir.display())),
    };
    for entry in entries {
        let path = entry.map_err(|e| e.to_string())?.path();
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if !(name.starts_with("unit-") && name.ends_with(".json")) {
            continue;
        }
        let Ok(doc) = read_json::<RelationUnitDocument>(&path) else {
            ignored += 1;
            continue;
        };
        let belongs = doc.schema_version == 1
            && doc.params_digest == digest
            && same_curve(doc.degree, doc.curve_a, doc.subfield, doc.curve_b, c)
            && doc.spec == *spec
            && doc.summands == summands
            && name == format!("unit-{:05}.json", doc.unit);
        if !belongs || units.contains_key(&doc.unit) {
            ignored += 1;
            continue;
        }
        units.insert(doc.unit, doc);
    }
    Ok((units, ignored))
}

/// Collect one work unit and write its file atomically.
#[allow(clippy::too_many_arguments)]
fn collect_unit(
    c: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    ic: &KoblitzIcOptions,
    pair: &mut Option<PairSumTable>,
    p: &WorkflowParams,
    digest: &str,
    spec: &FactorBaseSpec,
    rel_dir: &Path,
    unit: usize,
) -> Result<RelationUnitDocument, String> {
    if ic.strategy == DecompositionStrategy::PairTable && pair.is_none() {
        *pair = Some(PairSumTable::build(c, fb).ok_or("field too wide for the pair table")?);
    }
    let collector = RelationCollector::with_pair_table(c, fb, ic, pair.as_ref())
        .ok_or("factor base cannot decompose with this summand count")?;
    let work = RelationWorkUnit {
        seed: p.seed,
        start: (unit as u64)
            .checked_mul(p.collection.unit_trials)
            .ok_or("work unit range overflows")?,
        count: p.collection.unit_trials,
    };
    let (relations, report) = collector.collect(work);
    let doc = RelationUnitDocument {
        schema_version: 1,
        params_digest: digest.to_string(),
        degree: c.n,
        curve_a: c.a,
        subfield: c.k,
        curve_b: c.b_index,
        spec: spec.clone(),
        summands: p.summands,
        unit,
        seed: work.seed,
        start: work.start,
        count: work.count,
        relations,
        elapsed_seconds: report.elapsed_seconds,
    };
    write_atomic(&unit_path(rel_dir.parent().unwrap_or(rel_dir), unit), &doc)?;
    Ok(doc)
}
