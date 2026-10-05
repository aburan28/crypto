//! The method registry: every ECDLP algorithm `ecbench` can run, behind
//! one interface, charged in one unit.
//!
//! A method is an id from [`registry`] plus string parameters.  Defaults
//! are written into the parameters before the method id is hashed, so an
//! omitted default and an explicit one are the same method, and an
//! unknown parameter is refused rather than ignored.
//!
//! Every method runs on the instance's [`CountedGroup`], so an addition
//! costs the same in every row: rho's walk step, a BSGS giant step, a
//! kangaroo jump and an index-calculus relation trial are all charged
//! through `GroupOps`.  Native work the unit does not price is counted
//! under a name ending `_uncharged` and listed in
//! [`SolveReport::unpriced`]; a report with any unpriced work is a lower
//! bound and says so.

use std::collections::BTreeMap;
use std::time::Instant;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::{derive_u64, sha256_hex, short_id};
use crate::cryptanalysis::ecbench::generic::{
    bsgs_interleaved, bsgs_negation, bsgs_textbook, kangaroo, GenericOutcome,
};
use crate::cryptanalysis::ecbench::workload::{CurveFacts, Instance};
use crate::cryptanalysis::ecbench_large_prime::{
    self as large_prime, SolveConfig as LargePrimeConfig,
};
use crate::cryptanalysis::ic_boundary::{
    rho_cap, rho_reference, rho_walk_with, signed_frobenius_rho_tuned, BinaryGroup, BinaryInstance,
    Calibration, CountedGroup, GroupOps, NegationClasses, PhaseCost, PointClasses, PrimeInstance,
    PrimePoint, RhoResult, RhoWalk,
};
use crate::cryptanalysis::ic_framework::plugins::{
    BinarySubspaceBase, CompactOrbitScanBase, DescentAlgebraicOracle, FrobeniusMitmOracle,
    GlvOrbitBase, KoblitzOrbitBase, MitmOracle, PrimeAbscissaBase, SubtractOracle,
};
use crate::cryptanalysis::ic_framework::shared_rank::{
    run_shared_rank_targets, SharedRankSpec, SharedTargetSpec,
};
use crate::cryptanalysis::ic_framework::solvers::solver_by_name;
use crate::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, Targets,
};
use crate::cryptanalysis::ic_framework::{run_pipeline, PipelineSpec, RunReport};
use crate::cryptanalysis::ic_measurement::{self as measurement, Phase, Snapshot};
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_strong_rho::{
    RawPoint as StrongPoint, StrongRho, StrongRhoCharges, StrongRhoParams,
};
use rand::{rngs::StdRng, SeedableRng};

/// A parameter a method reads.  `default: None` means required: a value
/// that changes the search is never left to a producer default.
#[derive(Clone, Copy, Debug)]
pub struct ParamDecl {
    pub name: &'static str,
    pub default: Option<&'static str>,
    pub help: &'static str,
}

/// Which curves a method runs on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Applies {
    Any,
    KoblitzOnly,
    BinaryOnly,
}

/// One registered method.
#[derive(Clone, Copy, Debug)]
pub struct MethodDecl {
    pub id: &'static str,
    /// `rho`, `bsgs`, `kangaroo` or `ic`.
    pub family: &'static str,
    pub summary: &'static str,
    /// The code the method runs, for the record.
    pub entry: &'static str,
    pub applies: Applies,
    pub params: &'static [ParamDecl],
}

const RHO_PARAMS: &[ParamDecl] = &[ParamDecl {
    name: "cap_multiple",
    default: Some("64"),
    help: "step budget as a multiple of √r before the run counts as exhausted",
}];

/// Every method `ecbench` knows.  Adding one is an entry here and an arm
/// in [`solve`]; nothing else changes.
pub fn registry() -> &'static [MethodDecl] {
    &[
        MethodDecl {
            id: "rho.frozen_reference",
            family: "rho",
            summary: "the frozen r-adding walk every ledger row was priced against through §17 (A = 1); the historical \"before\" mark",
            entry: "ic_boundary::rho_reference",
            applies: Applies::Any,
            params: RHO_PARAMS,
        },
        MethodDecl {
            id: "rho.plain",
            family: "rho",
            summary: "tuned r-adding walk on points (A = 1), distinguished points, stride starts",
            entry: "ic_boundary::rho_walk_with(PointClasses, RhoWalk::plain())",
            applies: Applies::Any,
            params: RHO_PARAMS,
        },
        MethodDecl {
            id: "rho.negation",
            family: "rho",
            summary: "tuned walk on {P, −P} (A = 2) with Wiener–Zuccherato look-ahead and cycle escape",
            entry: "ic_boundary::rho_walk_with(NegationClasses, RhoWalk::negation())",
            applies: Applies::Any,
            params: RHO_PARAMS,
        },
        MethodDecl {
            id: "rho.signed_frobenius",
            family: "rho",
            summary: "tuned walk on signed Frobenius orbits {±φ^t P} (A = 2n), Koblitz curves only",
            entry: "ic_boundary::signed_frobenius_rho_tuned",
            applies: Applies::KoblitzOnly,
            params: RHO_PARAMS,
        },
        MethodDecl {
            id: "rho.signed_frobenius_strong",
            family: "rho",
            summary: "the strong single-target reference of the IC measurement rules (docs/ic/boundary_targets.json, since 2026-10-01): lockstep walks with batched inversion on signed Frobenius orbits (A = 2n), Koblitz curves only",
            entry: "koblitz_strong_rho::StrongRho::solve",
            applies: Applies::KoblitzOnly,
            params: &[
                ParamDecl {
                    name: "lanes",
                    default: Some("32"),
                    help: "walks advanced in lockstep (batch-inversion width); the reference measured 32",
                },
                ParamDecl {
                    name: "dp_bits",
                    default: Some("4"),
                    help: "distinguished-point bits; the reference measured 4",
                },
                ParamDecl {
                    name: "step_cap_factor",
                    default: Some("2000"),
                    help: "give up after this multiple of the ideal step count",
                },
            ],
        },
        MethodDecl {
            id: "bsgs.textbook",
            family: "bsgs",
            summary: "baby table of m = ⌈√r⌉ steps first, then giant steps; worst 2√r, mean 1.5√r",
            entry: "ecbench::generic::bsgs_textbook",
            applies: Applies::Any,
            params: &[],
        },
        MethodDecl {
            id: "bsgs.interleaved",
            family: "bsgs",
            summary: "Pollard's interleaved baby and giant steps, two tables; mean (4/3)√r",
            entry: "ecbench::generic::bsgs_interleaved",
            applies: Applies::Any,
            params: &[],
        },
        MethodDecl {
            id: "bsgs.negation",
            family: "bsgs",
            summary: "baby table keyed by {P, −P}, giant stride 2m + 1; mean about √r",
            entry: "ecbench::generic::bsgs_negation",
            applies: Applies::Any,
            params: &[],
        },
        MethodDecl {
            id: "kangaroo.vow",
            family: "kangaroo",
            summary: "van Oorschot–Wiener tame/wild kangaroo on [0, r), distinguished points; mean about 2√r",
            entry: "ecbench::generic::kangaroo",
            applies: Applies::Any,
            params: RHO_PARAMS,
        },
        MethodDecl {
            id: "ic.pipeline",
            family: "ic",
            summary: "index calculus through ic_framework::run_pipeline: factor base, oracle set-up, relations, linear algebra, verification, all charged",
            entry: "ic_framework::run_pipeline",
            applies: Applies::Any,
            params: &[
                ParamDecl {
                    name: "factor_base",
                    default: None,
                    help: "plug-in spec: prime-abscissa:size=N, glv-orbit:size=N, binary-subspace:dimension=D, koblitz-orbit:divisor=1;2, compact-orbit-scan:columns=N,raw_x_cap=M",
                },
                ParamDecl {
                    name: "oracle",
                    default: None,
                    help: "subtract, mitm, mitm-frobenius, mitm-frobenius-counted or descent-algebraic, with :m=2|3",
                },
                ParamDecl {
                    name: "solver",
                    default: Some(""),
                    help: "descent-algebraic only: buchberger-f2, sat-cdcl, exhaustive (with :k=v options)",
                },
                ParamDecl {
                    name: "linalg",
                    default: Some("incremental-gauss"),
                    help: "incremental-gauss (stop at target pin), incremental-gauss-full-rank (all columns independent), incremental-gauss-full-rank-checked (also point-check each base column), or structured-gauss",
                },
                ParamDecl {
                    name: "targets",
                    default: Some("walk"),
                    help: "random ([a]G+[b]Q per trial) or walk (one addition per trial)",
                },
                ParamDecl {
                    name: "max_trials",
                    default: Some("100000000"),
                    help: "relation trials before the run counts as exhausted",
                },
                ParamDecl {
                    name: "solver_budget_seconds",
                    default: Some("0"),
                    help: "per-call wall budget for an algebraic solver; nonzero makes the run nondeterministic",
                },
            ],
        },
        MethodDecl {
            id: "ic.large_prime",
            family: "ic",
            summary: "bounded exact m-summand binary index calculus with zero, single, or double large primes and exact modular elimination",
            entry: "ecbench_large_prime::solve",
            applies: Applies::BinaryOnly,
            params: &[
                ParamDecl {
                    name: "small_dimension",
                    default: None,
                    help: "dimension of the nested small x-coordinate subspace",
                },
                ParamDecl {
                    name: "envelope_dimension",
                    default: None,
                    help: "dimension of the full x-coordinate factor-base envelope",
                },
                ParamDecl {
                    name: "summands",
                    default: None,
                    help: "exact summand count, or n-1 (never silently downscaled)",
                },
                ParamDecl {
                    name: "large_primes",
                    default: None,
                    help: "maximum accepted large-prime columns: 0, 1, or 2",
                },
                ParamDecl {
                    name: "max_trials",
                    default: None,
                    help: "relation trials before an exhausted result",
                },
                ParamDecl {
                    name: "max_states",
                    default: None,
                    help: "hard cap on exact meet-in-the-middle combination states",
                },
            ],
        },
        MethodDecl {
            id: "ic.shared_rank",
            family: "ic",
            summary: "target-blind full-rank Koblitz relation log followed by point-only signed-Frobenius m3 target descent; native work remains visible and unpriced",
            entry: "ic_framework::shared_rank::run_shared_rank_targets",
            applies: Applies::KoblitzOnly,
            params: &[
                ParamDecl {
                    name: "factor_base",
                    default: None,
                    help: "compact-orbit-scan:columns=N,raw_x_cap=M; exact base construction",
                },
                ParamDecl {
                    name: "rank_seed",
                    default: None,
                    help: "seed of the target-blind relation search",
                },
                ParamDecl {
                    name: "rank_max_trials",
                    default: None,
                    help: "target-blind relation trial cap",
                },
                ParamDecl {
                    name: "target_max_attempts",
                    default: None,
                    help: "direct Q query followed by at most this many total residual attempts",
                },
            ],
        },
    ]
}

/// The leading-order expected `S` of a generic method on a curve whose
/// generic automorphism group has order `a_available`, or `None` when
/// the method has no such constant (index calculus).  Derived, not
/// measured: these are the boundaries a calibration run is read against.
///
/// - rho on points: `√(π/2)` (the birthday bound on `r` points);
///   on classes of size `A`: `√(π/2A)` (negation `A = 2`, signed
///   Frobenius `A = 2m`).
/// - BSGS, textbook: `m = ⌈√r⌉` baby steps, then a uniform target's
///   giant steps, mean `m/2`: `1.5`.  Interleaved: `E[2·max(i, j)]` for
///   uniform `i, j < m`, `(4/3)`.  Negation-folded: `m = √r/2` baby
///   steps and a mean of `√r/2` giant steps of stride `2m + 1`: `1`.
/// - Kangaroo (van Oorschot–Wiener, one tame and one wild, mean jump
///   `√r/2` on the interval `[0, r)`): `2`.
///
/// Every one omits `O(log r)` set-up, so measured `S` approaches it from
/// above as `r` grows.
pub fn expected_s(id: &str, a_available: u32) -> Option<f64> {
    use std::f64::consts::PI;
    match id {
        "rho.frozen_reference" | "rho.plain" => Some((PI / 2.0).sqrt()),
        "rho.negation" => Some((PI / 4.0).sqrt()),
        "rho.signed_frobenius" | "rho.signed_frobenius_strong" => {
            Some((PI / (2.0 * a_available.max(1) as f64)).sqrt())
        }
        "bsgs.textbook" => Some(1.5),
        "bsgs.interleaved" => Some(4.0 / 3.0),
        "bsgs.negation" => Some(1.0),
        "kangaroo.vow" => Some(2.0),
        _ => None,
    }
}

/// A method as a spec names it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MethodSpec {
    pub id: String,
    #[serde(default)]
    pub params: BTreeMap<String, String>,
}

/// A method with its defaults written out and its identity computed.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ResolvedMethod {
    pub id: String,
    pub family: String,
    pub params: BTreeMap<String, String>,
    /// `ECM1h` + 12 hex of `{id, params}`.  The configuration, not the
    /// code: the binary's hash is recorded beside it.
    pub method_id: String,
    pub method_sha256: String,
    pub entry: String,
}

pub fn decl(id: &str) -> Option<&'static MethodDecl> {
    registry().iter().find(|m| m.id == id)
}

/// Fill defaults, refuse unknown and missing parameters, hash.
pub fn resolve(spec: &MethodSpec) -> Result<ResolvedMethod, String> {
    let d = decl(&spec.id).ok_or_else(|| {
        let ids: Vec<&str> = registry().iter().map(|m| m.id).collect();
        format!("unknown method `{}`; known: {}", spec.id, ids.join(", "))
    })?;
    for k in spec.params.keys() {
        if !d.params.iter().any(|p| p.name == k) {
            return Err(format!("method `{}` has no parameter `{k}`", d.id));
        }
    }
    let mut params = BTreeMap::new();
    for p in d.params {
        match (spec.params.get(p.name), p.default) {
            (Some(v), _) => {
                params.insert(p.name.to_string(), v.clone());
            }
            (None, Some(def)) => {
                params.insert(p.name.to_string(), def.to_string());
            }
            (None, None) => {
                return Err(format!(
                    "method `{}` needs `{}` ({}); it has no default because it changes the search",
                    d.id, p.name, p.help
                ))
            }
        }
    }
    let (method_id, method_sha256) = short_id(
        "ECM1",
        &json!({"schema": "ecbench.method/v1", "id": d.id, "params": params}),
    )?;
    Ok(ResolvedMethod {
        id: d.id.into(),
        family: d.family.into(),
        params,
        method_id,
        method_sha256,
        entry: d.entry.into(),
    })
}

/// One phase of a solve, in the unit and natively.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PhaseRecord {
    pub name: String,
    pub adds: u64,
    pub doubles: u64,
    pub scalar_mults: u64,
    /// Group-addition equivalents charged to this phase.
    pub gae: f64,
    /// In-process wall time of the phase where the producer clocks it.
    pub wall_ns: Option<u64>,
    #[serde(default)]
    pub native: BTreeMap<String, u64>,
}

/// The factor base an index-calculus run used, identified by its points.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FactorBaseFacts {
    /// `FB1h` + 12 hex of the curve, family, parameters and the sorted
    /// point keys: two runs with the same id used the same points.
    pub fb_id: String,
    pub fb_sha256: String,
    pub family: String,
    pub params: BTreeMap<String, String>,
    pub description: String,
    pub signed_points: u64,
    pub abscissae: u64,
    pub columns: u64,
    pub dimension: Option<u32>,
    /// SHA-256 of the sorted group keys of every point.
    pub points_sha256: String,
}

/// The one-target online window: from the first target-dependent
/// computation to the verified recovery, with its exclusive phases, as the
/// IC measurement rules define it (AGENTS.md "IC measurements").
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OnlineWindow {
    pub wall_ns: u64,
    /// Exclusive phases inside the window, ns, under the claim schema's
    /// names: `target_query`, `target_PDP`, `target_relation_check`,
    /// `target_descent`, `target_recovery_check` for index calculus;
    /// `rho_solve` (and `recovery_check` where a walk marks it) for a
    /// generic method.  They sum to `wall_ns`.
    pub phases_ns: BTreeMap<String, u64>,
    /// Phases of the window this method never entered, each with the
    /// reason its cost is zero by construction.
    pub zero_phases: BTreeMap<String, String>,
    pub start_event: String,
    pub stop_event: String,
    pub included_stages: Vec<String>,
    /// How the method's work maps onto the stage names.
    pub mapping: String,
}

/// What an algebraic or SAT decomposition oracle's solver did over a
/// run: the shape of the system it solved, the work it counted in its
/// own unit, and the matrix and search statistics the ICMS registry
/// names (`docs/ic/measurement/registry.json`, `pdp_metrics`).
///
/// Informational: nothing here is priced into `S` beyond what the
/// relations phase already carries as `solver_*` counters, and the audit
/// replay does not compare it (a budgeted solver is nondeterministic by
/// construction, and the block records exactly that).  A table oracle
/// has no solver and the record's block is `null`.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SolverStats {
    pub name: String,
    pub calls: u64,
    pub ops: u64,
    pub op_unit: String,
    pub wall_ns: u64,
    pub peak_bytes: u64,
    /// Calls that ran out of budget, counted apart from refutations.
    pub budget_exceeded: u64,
    /// The system every call solved on this base: variables, equations,
    /// their degrees, and the semi-regular degree such a system reaches
    /// without exploitable structure.
    pub n_vars: Option<u64>,
    pub n_equations: Option<u64>,
    pub degrees: Vec<u32>,
    pub semi_regular_degree: Option<u32>,
    pub solving_degree_mean: Option<f64>,
    pub solving_degree_max: Option<u32>,
    /// Mean solving degree over the semi-regular bound.
    pub degree_over_bound: Option<f64>,
    /// Macaulay (or F4/F5) matrix figures, the maximum over calls where
    /// the solver reports a maximum and the sum where it reports a sum.
    pub macaulay_rows: Option<u64>,
    pub macaulay_columns: Option<u64>,
    pub macaulay_degree: Option<u32>,
    pub macaulay_rank: Option<u64>,
    /// SAT figures: variables and clauses are the maximum over calls,
    /// conflicts, decisions, propagations, restarts and learnt clauses
    /// are totals over the run.
    pub sat_variables: Option<u64>,
    pub sat_clauses: Option<u64>,
    pub sat_conflicts: Option<u64>,
    pub sat_decisions: Option<u64>,
    pub sat_propagations: Option<u64>,
    pub sat_restarts: Option<u64>,
    pub sat_learnt_clauses: Option<u64>,
    /// Every solver-specific counter, verbatim.
    pub extra: BTreeMap<String, u64>,
}

impl SolverStats {
    /// Build the block from the framework's run report.
    pub fn from_report(
        solver: &crate::cryptanalysis::ic_framework::SolverReport,
        system: Option<&crate::cryptanalysis::ic_framework::stages::SystemShape>,
    ) -> Self {
        let x = &solver.extra;
        let get = |keys: &[&str]| keys.iter().find_map(|k| x.get(*k).copied());
        let sat = solver.op_unit == "conflicts";
        SolverStats {
            name: solver.name.clone(),
            calls: solver.calls,
            ops: solver.ops,
            op_unit: solver.op_unit.clone(),
            wall_ns: solver.wall_ns,
            peak_bytes: solver.peak_bytes,
            budget_exceeded: solver.budget_exceeded,
            n_vars: system.map(|s| s.n_vars as u64),
            n_equations: system.map(|s| s.n_equations as u64),
            degrees: system.map(|s| s.degrees.clone()).unwrap_or_default(),
            semi_regular_degree: solver
                .semi_regular_degree
                .or_else(|| system.and_then(|s| s.semi_regular_degree)),
            solving_degree_mean: solver.solving_degree_mean,
            solving_degree_max: (solver.solving_degree_max > 0)
                .then_some(solver.solving_degree_max),
            degree_over_bound: solver.degree_over_bound,
            macaulay_rows: get(&["matrix_rows_max", "macaulay_rows"]),
            macaulay_columns: get(&["matrix_cols_max", "macaulay_cols"]),
            macaulay_degree: get(&["max_degree", "max_poly_degree"])
                .map(|d| d as u32)
                .or((solver.solving_degree_max > 0).then_some(solver.solving_degree_max)),
            macaulay_rank: get(&["matrix_rank_max", "macaulay_rank"]),
            sat_variables: if sat { get(&["variables_max"]) } else { None },
            sat_clauses: if sat { get(&["clauses_max"]) } else { None },
            sat_conflicts: sat.then_some(solver.ops),
            sat_decisions: if sat { get(&["decisions"]) } else { None },
            sat_propagations: if sat { get(&["propagations"]) } else { None },
            sat_restarts: if sat { get(&["restarts"]) } else { None },
            sat_learnt_clauses: if sat { get(&["learnt_clauses"]) } else { None },
            extra: x.clone(),
        }
    }
}

/// What a solve reports.  The runner adds verification and timing.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SolveReport {
    pub recovered: Option<u64>,
    /// The budget ran out before an answer.
    pub exhausted: bool,
    pub phases: Vec<PhaseRecord>,
    pub total_gae: f64,
    /// Automorphism group order the method's walk or table uses.
    pub automorphisms_used: u32,
    /// Counters by name; shape descriptors and `*_uncharged` work.
    pub counters: BTreeMap<String, u64>,
    /// Work counted but not priced in the unit.  Nonempty makes the
    /// total a lower bound.
    pub unpriced: Vec<String>,
    /// The same inputs give the same counts.  False when any price came
    /// from host wall time or a solver ran under a wall-clock budget.
    pub deterministic: bool,
    pub nondeterminism: Vec<String>,
    /// Wall time of the algorithm alone, inside the measured process.
    pub solve_wall_ns: u64,
    pub factor_base: Option<FactorBaseFacts>,
    /// Method-specific detail kept whole (the IC report's structure).
    pub detail: Value,
    /// The one-target online window, when the phase clock produced one.
    #[serde(default)]
    pub online: Option<OnlineWindow>,
    #[serde(default)]
    pub online_error: Option<String>,
    /// The decomposition solver's statistics, for algebraic and SAT
    /// oracles; `None` for table oracles and generic methods.
    #[serde(default)]
    pub solver: Option<SolverStats>,
}

fn ops_phase(name: &str, ops: GroupOps) -> PhaseRecord {
    PhaseRecord {
        name: name.into(),
        adds: ops.adds,
        doubles: ops.doubles,
        scalar_mults: ops.scalar_mults,
        gae: ops.gae(),
        wall_ns: None,
        native: BTreeMap::new(),
    }
}

fn unpriced_of(counters: &BTreeMap<String, u64>) -> Vec<String> {
    counters
        .iter()
        .filter(|(k, v)| k.ends_with("_uncharged") && **v > 0)
        .map(|(k, _)| k.clone())
        .collect()
}

fn param_u64(m: &ResolvedMethod, name: &str) -> Result<u64, String> {
    m.params[name]
        .parse()
        .map_err(|_| format!("parameter `{name}` is not an integer: `{}`", m.params[name]))
}

fn generic_report(out: GenericOutcome, wall: u64, a: u32) -> SolveReport {
    let total = out.setup.gae() + out.search.gae();
    SolveReport {
        recovered: out.recovered,
        exhausted: out.exhausted,
        phases: vec![
            ops_phase("setup", out.setup),
            ops_phase("search", out.search),
        ],
        total_gae: total,
        automorphisms_used: a,
        unpriced: unpriced_of(&out.counters),
        counters: out.counters,
        deterministic: true,
        nondeterminism: vec![],
        solve_wall_ns: wall,
        factor_base: None,
        detail: Value::Null,
        online: None,
        online_error: None,
        solver: None,
    }
}

/// Split a tuned rho run into set-up, walk and verification from its own
/// counters; the frozen walk has none and stays one phase.
fn rho_report(res: RhoResult, wall: u64) -> SolveReport {
    let c = &res.counters;
    let get = |k: &str| c.get(k).copied().unwrap_or(0);
    let phases = if c.is_empty() {
        vec![ops_phase("search", res.group_ops)]
    } else {
        let setup = get("setup_additions") + get("setup_doublings");
        let verify = get("verification_operations");
        let total = res.group_ops.adds + res.group_ops.doubles;
        let search = total.saturating_sub(setup + verify);
        let mk = |name: &str, gae: u64| PhaseRecord {
            name: name.into(),
            gae: gae as f64,
            ..Default::default()
        };
        let mut s = mk("setup", setup);
        s.adds = get("setup_additions");
        s.doubles = get("setup_doublings");
        let mut w = mk("search", search);
        w.scalar_mults = res.group_ops.scalar_mults;
        vec![s, w, mk("internal_verification", verify)]
    };
    let mut counters = res.counters.clone();
    counters.insert("steps".into(), res.steps);
    counters.insert("walks".into(), res.walks);
    counters.insert("distinguished_points".into(), res.distinguished_points);
    SolveReport {
        recovered: res.recovered,
        exhausted: res.recovered.is_none(),
        total_gae: res.gae,
        automorphisms_used: res.automorphisms,
        unpriced: unpriced_of(&counters),
        counters,
        phases,
        deterministic: true,
        nondeterminism: vec![],
        solve_wall_ns: wall,
        factor_base: None,
        detail: json!({"method": res.method, "expected_steps": res.expected_steps}),
        online: None,
        online_error: None,
        solver: None,
    }
}

/// Run the generic and rho methods on any counted group.
fn solve_generic<G: CountedGroup>(
    m: &ResolvedMethod,
    g: &G,
    gen: G::Elt,
    target: G::Elt,
    r: u64,
    seed: u64,
) -> Result<SolveReport, String> {
    let cap = || -> Result<u64, String> { Ok(rho_cap(r, param_u64(m, "cap_multiple")? as f64)) };
    let t = Instant::now();
    // Everything these methods do depends on the target (rho's jump table
    // is [a]G + [b]Q), so the whole solve is the online window.
    measurement::begin_online(Phase::RhoSolve);
    let rep = match m.id.as_str() {
        "rho.frozen_reference" => {
            let res = rho_reference(g, gen, target, r, seed, cap()?);
            rho_report(res, t.elapsed().as_nanos() as u64)
        }
        "rho.plain" => {
            let res = rho_walk_with(
                g,
                &PointClasses,
                gen,
                target,
                r,
                seed,
                cap()?,
                RhoWalk::plain(),
            );
            rho_report(res, t.elapsed().as_nanos() as u64)
        }
        "rho.negation" => {
            let res = rho_walk_with(
                g,
                &NegationClasses { r },
                gen,
                target,
                r,
                seed,
                cap()?,
                RhoWalk::negation(),
            );
            rho_report(res, t.elapsed().as_nanos() as u64)
        }
        "bsgs.textbook" => {
            let o = bsgs_textbook(g, gen, target, r);
            generic_report(o, t.elapsed().as_nanos() as u64, 1)
        }
        "bsgs.interleaved" => {
            let o = bsgs_interleaved(g, gen, target, r);
            generic_report(o, t.elapsed().as_nanos() as u64, 1)
        }
        "bsgs.negation" => {
            let o = bsgs_negation(g, gen, target, r);
            generic_report(o, t.elapsed().as_nanos() as u64, 2)
        }
        "kangaroo.vow" => {
            let o = kangaroo(g, gen, target, r, seed, cap()?);
            generic_report(o, t.elapsed().as_nanos() as u64, 1)
        }
        other => return Err(format!("`{other}` is not a generic method")),
    };
    measurement::end_online();
    Ok(rep)
}

fn parse_hex(s: &str) -> Result<u64, String> {
    u64::from_str_radix(s.trim_start_matches("0x"), 16).map_err(|_| format!("bad hex `{s}`"))
}

/// Solve for `target` with method `m` on `inst`, under algorithm seed
/// `seed`.  The planted logarithm is not an argument.
pub fn solve(
    m: &ResolvedMethod,
    inst: &Instance,
    curve: &CurveFacts,
    target: &[String; 2],
    seed: u64,
) -> Result<SolveReport, String> {
    let d = decl(&m.id).ok_or("unregistered method")?;
    if d.applies == Applies::KoblitzOnly && curve.family != "koblitz" {
        return Err(format!(
            "`{}` runs on Koblitz curves only; {} is {}",
            m.id, curve.slug, curve.family
        ));
    }
    if d.applies == Applies::BinaryOnly && !matches!(inst, Instance::Binary(_)) {
        return Err(format!(
            "`{}` runs on binary curves only; {} is {}",
            m.id, curve.slug, curve.family
        ));
    }
    let (tx, ty) = (parse_hex(&target[0])?, parse_hex(&target[1])?);
    // The phase clock: the methods open and close the online window, the
    // session turns it into exclusive phase times.
    let session = measurement::Session::begin().ok();
    let mut rep = solve_inner(m, inst, curve, d, tx, ty, seed)?;
    match session.map(|s| s.finish()) {
        Some(Ok(snap)) if rep.online.is_none() => {
            rep.online = online_window(d.family, &m.id, &snap)
        }
        Some(Ok(_)) => {} // A method with its own five-phase online clock supplied the window.
        Some(Err(e)) => rep.online_error = Some(e.to_string()),
        None => rep.online_error = Some("a measurement session was already active".into()),
    }
    Ok(rep)
}

fn solve_inner(
    m: &ResolvedMethod,
    inst: &Instance,
    curve: &CurveFacts,
    d: &MethodDecl,
    tx: u64,
    ty: u64,
    seed: u64,
) -> Result<SolveReport, String> {
    match (inst, d.family) {
        (Instance::Binary(i), _) if m.id == "ic.large_prime" => {
            solve_large_prime(m, i, curve, FastPoint::affine(tx, ty), seed)
        }
        (Instance::Binary(i), _) if m.id == "ic.shared_rank" => {
            solve_shared_rank(m, i, curve, FastPoint::affine(tx, ty), seed)
        }
        (Instance::Prime(i), "ic") => solve_ic_prime(m, i, curve, PrimePoint::affine(tx, ty), seed),
        (Instance::Binary(i), "ic") => {
            solve_ic_binary(m, i, curve, FastPoint::affine(tx, ty), seed)
        }
        (Instance::Binary(i), _) if m.id == "rho.signed_frobenius_strong" => {
            solve_strong(m, i, FastPoint::affine(tx, ty), seed)
        }
        (Instance::Binary(i), _) if m.id == "rho.signed_frobenius" => {
            let t = Instant::now();
            measurement::begin_online(Phase::RhoSolve);
            let res = signed_frobenius_rho_tuned(
                i,
                FastPoint::affine(tx, ty),
                seed,
                rho_cap(i.r, param_u64(m, "cap_multiple")? as f64),
            )
            .ok_or("not a Koblitz instance")?;
            measurement::end_online();
            Ok(rho_report(res, t.elapsed().as_nanos() as u64))
        }
        (Instance::Prime(i), _) => solve_generic(
            m,
            &i.curve,
            i.generator_point(),
            PrimePoint::affine(tx, ty),
            i.r,
            seed,
        ),
        (Instance::Binary(i), _) => {
            let g = BinaryGroup(&i.fast);
            solve_generic(m, &g, i.generator, FastPoint::affine(tx, ty), i.r, seed)
        }
    }
}

fn solve_large_prime(
    m: &ResolvedMethod,
    inst: &BinaryInstance,
    curve: &CurveFacts,
    target: FastPoint,
    seed: u64,
) -> Result<SolveReport, String> {
    let summands = match m.params["summands"].as_str() {
        "n-1" => inst.n.checked_sub(1).ok_or("field degree has no n-1")?,
        value => value
            .parse::<u32>()
            .map_err(|_| format!("parameter `summands` is not an integer or n-1: `{value}`"))?,
    };
    let config = LargePrimeConfig {
        small_dimension: u32::try_from(param_u64(m, "small_dimension")?)
            .map_err(|_| "small_dimension does not fit u32")?,
        envelope_dimension: u32::try_from(param_u64(m, "envelope_dimension")?)
            .map_err(|_| "envelope_dimension does not fit u32")?,
        summands,
        max_large_primes: u8::try_from(param_u64(m, "large_primes")?)
            .map_err(|_| "large_primes does not fit u8")?,
        max_trials: param_u64(m, "max_trials")?,
        max_states: param_u64(m, "max_states")?,
        seed,
    };
    // The method returns an algebraic candidate; ecbench's runner verifies it.
    let report = large_prime::solve_candidate(inst, target, config)?;
    let mut phases = Vec::new();
    for name in ["factor_base", "oracle_setup", "relations", "verification"] {
        let mut phase = ops_phase(
            name,
            report
                .phase_group_ops
                .get(name)
                .copied()
                .unwrap_or_default(),
        );
        match name {
            "oracle_setup" => {
                phase.native.insert(
                    "combination_states_uncharged".into(),
                    report.counters.combination_states_uncharged,
                );
            }
            "relations" => {
                phase.native.insert(
                    "mitm_lookups_uncharged".into(),
                    report.counters.mitm_lookups_uncharged,
                );
                phase.native.insert(
                    "lp_merge_ops_uncharged".into(),
                    report.counters.lp_merge_ops_uncharged,
                );
                phase.native.insert(
                    "row_ops_uncharged".into(),
                    report.counters.row_ops_uncharged,
                );
            }
            _ => {}
        }
        phases.push(phase);
    }
    let total_gae = phases.iter().map(|p| p.gae).sum();
    let mut counters = BTreeMap::new();
    for (name, value) in [
        ("trials", report.counters.trials),
        ("decompositions", report.counters.decompositions),
        ("accepted_relations", report.counters.accepted_relations),
        (
            "eliminated_full_relations",
            report.counters.eliminated_full_relations,
        ),
        ("matrix_rank", report.counters.matrix_rank),
        ("large_prime_pivots", report.counters.lp_pivots),
        (
            "combination_states_uncharged",
            report.counters.combination_states_uncharged,
        ),
        (
            "mitm_lookups_uncharged",
            report.counters.mitm_lookups_uncharged,
        ),
        (
            "lp_merge_ops_uncharged",
            report.counters.lp_merge_ops_uncharged,
        ),
        ("row_ops_uncharged", report.counters.row_ops_uncharged),
    ] {
        counters.insert(name.into(), value);
    }
    for (i, value) in report.counters.relation_histogram.iter().enumerate() {
        counters.insert(format!("relations_with_{}_large_primes", i.min(3)), *value);
    }
    let mut point_bytes = Vec::with_capacity(report.factor_base.column_keys.len() * 8);
    for key in &report.factor_base.column_keys {
        point_bytes.extend_from_slice(&key.to_be_bytes());
    }
    let points_sha256 = sha256_hex(&point_bytes);
    let (fb_id, fb_sha256) = short_id(
        "FB1",
        &json!({
            "schema": "ecbench.factor_base/v1",
            "curve": curve.slug,
            "family": "binary-nested-large-prime",
            "params": m.params,
            "columns": report.factor_base.columns,
            "signed_points": report.factor_base.raw_points,
            "points_sha256": points_sha256,
        }),
    )?;
    let factor_base = FactorBaseFacts {
        fb_id,
        fb_sha256,
        family: "binary-nested-large-prime".into(),
        params: m.params.clone(),
        description: format!(
            "x < 2^{} nested at x < 2^{}, projected by cofactor and folded by negation",
            report.config.envelope_dimension, report.config.small_dimension
        ),
        signed_points: report.factor_base.raw_points as u64,
        abscissae: 1u64 << report.config.envelope_dimension,
        columns: report.factor_base.columns as u64,
        dimension: Some(report.config.envelope_dimension),
        points_sha256,
    };
    let detail =
        serde_json::to_value(&report).map_err(|e| format!("serialise large-prime report: {e}"))?;
    Ok(SolveReport {
        recovered: report.recovered,
        exhausted: report.exhausted,
        phases,
        total_gae,
        automorphisms_used: 2,
        counters,
        unpriced: vec![
            "combination_states_uncharged".into(),
            "mitm_lookups_uncharged".into(),
            "lp_merge_ops_uncharged".into(),
            "row_ops_uncharged".into(),
        ],
        deterministic: true,
        nondeterminism: vec![],
        solve_wall_ns: report.solve_wall_ns,
        factor_base: Some(factor_base),
        detail,
        online: None,
        online_error: None,
        solver: None,
    })
}

/// The online window from the phase clock's snapshot, mapped to the
/// claim schema's stage names.
fn online_window(family: &str, id: &str, snap: &Snapshot) -> Option<OnlineWindow> {
    let wall_ns = snap.online_wall_ns?;
    let get = |name: &str| snap.online_phases_ns.get(name).copied().flatten();
    let mut w = OnlineWindow {
        wall_ns,
        ..Default::default()
    };
    if family == "ic" {
        w.start_event = "first target-dependent query ([a]G + [b]Q, or the walk's first jump), after the factor base and the oracle's tables are built".into();
        w.stop_event = "[d]G = Q verified inside the pipeline".into();
        w.mapping = "classic index calculus: every relation carries the target, so query generation is target_query, decomposition attempts are target_PDP, the oracle's own witness checks are target_relation_check, the elimination that pins the target's column is target_descent, and the final [d]G = Q is target_recovery_check".into();
        for (clock, claim, why) in [
            ("target_query", "target_query", "no query was drawn"),
            ("target_pdp", "target_PDP", "no decomposition was attempted"),
            (
                "target_relation_check",
                "target_relation_check",
                "the oracle returns exact decompositions (pair-table hits on exact point keys, or solutions it checks under its own relation_check scope); this pipeline has no separate relation check, so any checking cost is inside target_PDP",
            ),
            ("target_descent", "target_descent", "no relation reached the elimination"),
            ("recovery_check", "target_recovery_check", "no candidate logarithm was pinned"),
        ] {
            match get(clock) {
                Some(ns) => {
                    w.phases_ns.insert(claim.into(), ns);
                }
                None => {
                    w.phases_ns.insert(claim.into(), 0);
                    w.zero_phases.insert(claim.into(), why.into());
                }
            }
            w.included_stages.push(claim.into());
        }
    } else {
        for (name, ns) in &snap.online_phases_ns {
            if let Some(ns) = ns {
                w.phases_ns.insert((*name).to_string(), *ns);
            }
        }
        let strong = id == "rho.signed_frobenius_strong";
        w.start_event = if strong {
            "first walk start [c]G + Q (the jump table [a]G is target-independent set-up)".into()
        } else {
            "first target-dependent operation (the jump table [a]G + [b]Q, or the first table step)"
                .into()
        };
        w.stop_event = "[d]G = Q verified inside the method".into();
        w.included_stages = match family {
            "rho" => vec!["walk".into(), "collision".into(), "recovery_check".into()],
            "kangaroo" => vec!["jumps".into(), "collision".into(), "recovery_check".into()],
            _ => vec![
                "baby_steps".into(),
                "giant_steps".into(),
                "recovery_check".into(),
            ],
        };
        w.mapping = format!("{family}: the whole target-dependent solve, one exclusive phase");
    }
    Some(w)
}

/// The strong reference: `koblitz_strong_rho` as the admissible fixture
/// drives it (`examples/koblitz_rho_fixture.rs`, backend `strong`), seeded
/// from the algorithm seed.  Group additions are exact; each scalar
/// multiplication is charged `1.5·log₂ r` additions, the convention
/// `ic_boundary::signed_frobenius_rho` uses.
fn solve_strong(
    m: &ResolvedMethod,
    inst: &BinaryInstance,
    target: FastPoint,
    seed: u64,
) -> Result<SolveReport, String> {
    let profile_online = std::env::var("ECBENCH_CALLGRIND_TARGET").as_deref() == Ok("1");
    let kc = inst.koblitz.as_ref().ok_or("not a Koblitz instance")?;
    let params = StrongRhoParams {
        lanes: param_u64(m, "lanes")?.max(1) as usize,
        dp_bits: param_u64(m, "dp_bits")? as u32,
        step_cap_factor: param_u64(m, "step_cap_factor")?,
    };
    if params.dp_bits >= 32 {
        return Err("dp_bits must be below 32".into());
    }
    let mul_gae = 1.5 * (inst.r as f64).log2();
    let t = Instant::now();
    let rho = StrongRho::new(kc);
    let mut charges = StrongRhoCharges::default();
    let jump_seed = derive_u64("ecbench.strong_rho.jumps", &[seed]);
    let jumps = rho.jumps(jump_seed, &mut charges);
    let setup_mults = charges.scalar_multiplications;
    let mut start_rng = StdRng::seed_from_u64(seed);
    let q = StrongPoint::from_binary(&inst.fast.lower(target));
    if profile_online {
        measurement::callgrind_dump(b"ecbench_before_online\0");
    }
    measurement::begin_online(Phase::RhoSolve);
    let outcome = rho.solve(q, &jumps, &mut start_rng, &params, charges);
    measurement::end_online();
    if profile_online {
        measurement::callgrind_dump(b"ecbench_online\0");
    }
    let wall = t.elapsed().as_nanos() as u64;
    let detail = json!({
        "reference": "koblitz_strong_rho::StrongRho, the single-target reference of docs/ic/boundary_targets.json since 2026-10-01",
        "lanes": params.lanes,
        "dp_bits": params.dp_bits,
        "jump_seed": jump_seed.to_string(),
        "scalar_multiplication_charge": "1.5·log2(r) additions each (ic_boundary::signed_frobenius_rho's convention)",
    });
    let Some(o) = outcome else {
        return Ok(SolveReport {
            recovered: None,
            exhausted: true,
            phases: vec![],
            total_gae: 0.0,
            automorphisms_used: 2 * inst.n,
            counters: BTreeMap::new(),
            unpriced: vec!["counts_lost_at_the_step_cap_uncharged".into()],
            deterministic: true,
            nondeterminism: vec![],
            solve_wall_ns: wall,
            factor_base: None,
            detail,
            online: None,
            online_error: None,
            solver: None,
        });
    };
    let c = o.charges;
    let search_mults = c.scalar_multiplications - setup_mults;
    let setup = PhaseRecord {
        name: "setup".into(),
        scalar_mults: setup_mults,
        gae: setup_mults as f64 * mul_gae,
        ..Default::default()
    };
    let search = PhaseRecord {
        name: "search".into(),
        adds: c.group_additions,
        scalar_mults: search_mults,
        gae: c.group_additions as f64 + search_mults as f64 * mul_gae,
        ..Default::default()
    };
    let mut counters = BTreeMap::new();
    for (k, v) in [
        ("walk_steps", o.walk_steps),
        ("walks", o.walks),
        ("fruitless_walks", o.fruitless),
        ("capped_walks", o.capped),
        ("wasted_merges", o.wasted_merges),
        ("table_entries", o.table_entries as u64),
        ("failed_collisions", c.failed_collisions),
        ("scalar_multiplications", c.scalar_multiplications),
        ("canonicalisations_uncharged", c.canonicalizations),
        ("partition_hashes_uncharged", c.partition_hashes),
        ("table_queries_uncharged", c.table_queries),
        ("table_inserts_uncharged", c.table_inserts),
    ] {
        counters.insert(k.to_string(), v);
    }
    Ok(SolveReport {
        recovered: Some(o.scalar),
        exhausted: false,
        total_gae: setup.gae + search.gae,
        phases: vec![setup, search],
        automorphisms_used: o.automorphisms as u32,
        unpriced: unpriced_of(&counters),
        counters,
        deterministic: true,
        nondeterminism: vec![],
        solve_wall_ns: wall,
        factor_base: None,
        detail,
        online: None,
        online_error: None,
        solver: None,
    })
}

// ── Index calculus ─────────────────────────────────────────────────

/// A calibration that prices only what the repository has pinned for
/// this curve, and nothing from the host: `ns_per_add = 1`, and every
/// unit without a pinned ratio left at zero (and reported unpriced).
/// The total is then a function of the counts alone.
fn pinned_only_calibration(regime: &str, slug: &str) -> Calibration {
    let mut c = Calibration {
        ns_per_add: 1.0,
        ns_per_double: 1.0,
        ..Default::default()
    };
    c.pin(regime, slug);
    c
}

/// Native counters the calibration prices, with the unit field that
/// prices each (`price_phase`'s table).
const PRICED_NATIVE: &[(&str, &str)] = &[
    ("sqrt_solves", "ns_per_sqrt"),
    ("as_solves", "ns_per_as_solve"),
    ("s4_pairs", "ns_per_s4_pair"),
    ("lookups", "ns_per_lookup"),
    ("row_ops", "ns_per_row_op"),
    ("legendre_symbols", "ns_per_legendre"),
    ("inversions", "ns_per_inversion"),
    ("frobenius_maps", "ns_per_frobenius"),
    ("canonicalisations", "ns_per_canon"),
    ("target_guard_probes", "ns_per_lookup"),
];

fn phase_of(name: &str, p: &PhaseCost) -> PhaseRecord {
    PhaseRecord {
        name: name.into(),
        adds: p.group_ops.adds,
        doubles: p.group_ops.doubles,
        scalar_mults: p.group_ops.scalar_mults,
        gae: p.gae,
        wall_ns: Some(p.wall_ns),
        native: p.native.clone(),
    }
}

fn ic_report(
    m: &ResolvedMethod,
    rep: RunReport,
    calib: &Calibration,
    fb: FactorBaseFacts,
    wall: u64,
    automorphisms: u32,
) -> SolveReport {
    let mut phases = vec![
        phase_of("factor_base", &rep.factor_base.cost),
        phase_of("oracle_setup", &rep.decomposition.setup),
        phase_of("relations", &rep.decomposition.cost),
        phase_of("linear_algebra", &rep.linear_algebra.cost),
        phase_of("verify", &rep.verify),
    ];
    let mut unpriced = Vec::new();
    let mut nondeterminism = Vec::new();
    // Solver work: run_pipeline prices it from this host's wall time,
    // which is not a count.  Take it back out and list it unpriced.
    let mut total = rep.total_gae;
    if let Some(s) = &rep.decomposition.solver {
        if s.gae > 0.0 && s.priced_by != "pinned" {
            phases[2].gae -= s.gae;
            total -= s.gae;
            unpriced.push(format!("solver_{}_uncharged", s.op_unit.replace(' ', "_")));
        }
        if param_u64(m, "solver_budget_seconds").unwrap_or(0) > 0 {
            nondeterminism.push("solver ran under a wall-clock budget".into());
        }
    }
    // Native counters with no pinned ratio for this curve.
    for (counter, unit) in PRICED_NATIVE {
        let present = phases
            .iter()
            .any(|p| p.native.get(*counter).copied().unwrap_or(0) > 0);
        if present && !calib.is_pinned(unit) {
            let name = format!("{counter}_uncharged");
            if !unpriced.contains(&name) {
                unpriced.push(name);
            }
        }
    }
    for phase in &phases {
        for (name, &count) in &phase.native {
            if count > 0 && name.ends_with("_uncharged") && !unpriced.contains(name) {
                unpriced.push(name.clone());
            }
        }
    }
    let mut counters = BTreeMap::new();
    counters.insert("targets_tried".into(), rep.decomposition.targets_tried);
    counters.insert("relations_found".into(), rep.decomposition.relations_found);
    counters.insert("matrix_rows".into(), rep.linear_algebra.rows);
    counters.insert("matrix_rank".into(), rep.linear_algebra.rank);
    counters.insert("matrix_dependent".into(), rep.linear_algebra.dependent);
    counters.insert("matrix_work".into(), rep.linear_algebra.work);
    let solver = rep
        .decomposition
        .solver
        .as_ref()
        .map(|sr| SolverStats::from_report(sr, rep.decomposition.system.as_ref()));
    let detail = json!({
        "label": rep.label,
        "decomposition": {
            "name": rep.decomposition.name,
            "summands": rep.decomposition.summands,
            "hit_rate": format!("{:.6e}", rep.decomposition.hit_rate),
            "system": rep.decomposition.system,
            "solver": rep.decomposition.solver,
        },
        "linear_algebra": {"name": rep.linear_algebra.name, "work_unit": rep.linear_algebra.work_unit},
        "calibration_pinned_units": calib.pinned_units,
        "framework_total_gae_before_solver_removal": rep.total_gae,
    });
    SolveReport {
        recovered: rep.recovered,
        exhausted: rep.exhausted,
        phases,
        total_gae: total,
        automorphisms_used: automorphisms,
        counters,
        unpriced,
        deterministic: nondeterminism.is_empty(),
        nondeterminism,
        solve_wall_ns: wall,
        factor_base: Some(fb),
        detail,
        online: None,
        online_error: None,
        solver,
    }
}

fn factor_base_facts<E: Copy>(
    slug: &str,
    spec: &str,
    params: &Params,
    fb: &crate::cryptanalysis::ic_boundary::FactorBase<E>,
    key: impl Fn(&E) -> u64,
) -> Result<FactorBaseFacts, String> {
    let mut keys: Vec<u64> = fb.points.iter().map(key).collect();
    keys.sort_unstable();
    let mut bytes = Vec::with_capacity(keys.len() * 8);
    for k in &keys {
        bytes.extend_from_slice(&k.to_be_bytes());
    }
    let points_sha256 = sha256_hex(&bytes);
    let family = spec.split(':').next().unwrap_or(spec).to_string();
    let (fb_id, fb_sha256) = short_id(
        "FB1",
        &json!({
            "schema": "ecbench.factor_base/v1",
            "curve": slug,
            "family": family,
            "params": params.0,
            "columns": fb.columns,
            "signed_points": fb.points.len(),
            "points_sha256": points_sha256,
        }),
    )?;
    Ok(FactorBaseFacts {
        fb_id,
        fb_sha256,
        family,
        params: params.0.clone(),
        description: fb.description.clone(),
        signed_points: fb.points.len() as u64,
        abscissae: fb.abscissae as u64,
        columns: fb.columns as u64,
        dimension: fb.dimension,
        points_sha256,
    })
}

fn pipeline_spec(m: &ResolvedMethod, seed: u64) -> Result<(PipelineSpec, String, String), String> {
    let (fb_name, fb_params) = Params::parse_spec(&m.params["factor_base"])?;
    let (or_name, or_params) = Params::parse_spec(&m.params["oracle"])?;
    let solver = m.params["solver"].clone();
    let (solver_name, solver_params) = if solver.is_empty() {
        (None, Params::default())
    } else {
        let (n, p) = Params::parse_spec(&solver)?;
        (Some(n), p)
    };
    let spec = PipelineSpec {
        factor_base: fb_name.clone(),
        factor_base_params: fb_params,
        oracle: or_name.clone(),
        oracle_params: or_params,
        solver: solver_name,
        solver_params,
        targets: Targets::parse(&m.params["targets"])?,
        linalg: m.params["linalg"].clone(),
        max_trials: param_u64(m, "max_trials")?,
        seed,
        solver_budget_seconds: param_u64(m, "solver_budget_seconds")?,
    };
    Ok((spec, fb_name, or_name))
}

fn solve_ic_prime(
    m: &ResolvedMethod,
    inst: &PrimeInstance,
    curve: &CurveFacts,
    target: PrimePoint,
    seed: u64,
) -> Result<SolveReport, String> {
    let (spec, fb_name, or_name) = pipeline_spec(m, seed)?;
    let abscissa = PrimeAbscissaBase { instance: inst };
    let orbit = GlvOrbitBase { instance: inst };
    let base: &dyn FactorBaseBuilder<_> = match fb_name.as_str() {
        "prime-abscissa" => &abscissa,
        "glv-orbit" => &orbit,
        other => {
            return Err(format!(
                "factor base `{other}` does not run on a prime curve"
            ))
        }
    };
    let ms = spec.oracle_params.u64_or("m", 2)? as u32;
    let mut subtract = SubtractOracle;
    let mut mitm = MitmOracle::new(ms);
    let oracle: &mut dyn DecompositionOracle<_> = match or_name.as_str() {
        "subtract" => &mut subtract,
        "mitm" => &mut mitm,
        other => return Err(format!("oracle `{other}` does not run on a prime curve")),
    };
    let ctx = InstanceCtx {
        group: &inst.curve,
        generator: inst.generator_point(),
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: curve.slug.clone(),
        field_degree: None,
    };
    let calib = pinned_only_calibration("prime", &curve.slug);
    let t = Instant::now();
    // The runner, not the pipeline, checks the answer: pass a planted
    // value that cannot match so `verified` carries no information.
    let rep = run_pipeline(&ctx, &spec, base, oracle, u64::MAX, &calib, None)?;
    let wall = t.elapsed().as_nanos() as u64;
    let mut scratch = GroupOps::default();
    let fb = base.build(&ctx, &spec.factor_base_params, &mut scratch)?;
    let facts = factor_base_facts(
        &curve.slug,
        &m.params["factor_base"],
        &spec.factor_base_params,
        &fb,
        |p| inst.curve.key(p),
    )?;
    Ok(ic_report(m, rep, &calib, facts, wall, 2))
}

fn solve_ic_binary(
    m: &ResolvedMethod,
    inst: &BinaryInstance,
    curve: &CurveFacts,
    target: FastPoint,
    seed: u64,
) -> Result<SolveReport, String> {
    let (spec, fb_name, or_name) = pipeline_spec(m, seed)?;
    let g = BinaryGroup(&inst.fast);
    let subspace = BinarySubspaceBase { instance: inst };
    let orbit = KoblitzOrbitBase { instance: inst };
    let compact = CompactOrbitScanBase { instance: inst };
    let base: &dyn FactorBaseBuilder<BinaryGroup> = match fb_name.as_str() {
        "binary-subspace" => &subspace,
        "koblitz-orbit" => &orbit,
        "compact-orbit-scan" => &compact,
        other => {
            return Err(format!(
                "factor base `{other}` does not run on a binary curve"
            ))
        }
    };
    let ms = spec.oracle_params.u64_or("m", 2)? as u32;
    let mut subtract = SubtractOracle;
    let mut mitm = MitmOracle::new(ms);
    let mut frob = FrobeniusMitmOracle::new(ms, inst);
    let mut frob_counted = FrobeniusMitmOracle::new_counted(ms, inst);
    let mut algebraic = match or_name.as_str() {
        "descent-algebraic" => {
            let name = spec
                .solver
                .as_deref()
                .ok_or("descent-algebraic needs `solver`")?;
            let budget = (spec.solver_budget_seconds > 0)
                .then(|| std::time::Duration::from_secs(spec.solver_budget_seconds));
            Some(DescentAlgebraicOracle::new(
                ms,
                inst,
                solver_by_name(name)?,
                spec.solver_params.clone(),
                budget,
            ))
        }
        _ => None,
    };
    let oracle: &mut dyn DecompositionOracle<BinaryGroup> = match or_name.as_str() {
        "subtract" => &mut subtract,
        "mitm" => &mut mitm,
        "mitm-frobenius" => &mut frob,
        "mitm-frobenius-counted" => &mut frob_counted,
        "descent-algebraic" => algebraic.as_mut().expect("built above"),
        other => return Err(format!("oracle `{other}` does not run on a binary curve")),
    };
    let ctx = InstanceCtx {
        group: &g,
        generator: inst.generator,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: curve.slug.clone(),
        field_degree: Some(inst.n),
    };
    let regime = if inst.koblitz.is_some() {
        "koblitz"
    } else {
        "char2"
    };
    let calib = pinned_only_calibration(regime, &curve.slug);
    let t = Instant::now();
    let rep = run_pipeline(&ctx, &spec, base, oracle, u64::MAX, &calib, None)?;
    let wall = t.elapsed().as_nanos() as u64;
    let mut scratch = GroupOps::default();
    let fb = base.build(&ctx, &spec.factor_base_params, &mut scratch)?;
    let facts = factor_base_facts(
        &curve.slug,
        &m.params["factor_base"],
        &spec.factor_base_params,
        &fb,
        |p| g.key(p),
    )?;
    let a = if inst.koblitz.is_some() {
        2 * inst.n
    } else {
        2
    };
    Ok(ic_report(m, rep, &calib, facts, wall, a))
}

/// The reusable, target-blind rank variant.  The library owns its exclusive
/// five-phase target clock; `ecbench` retains that clock instead of trying to
/// infer the online interval from the whole (mostly setup) process.
fn solve_shared_rank(
    m: &ResolvedMethod,
    inst: &BinaryInstance,
    curve: &CurveFacts,
    target: FastPoint,
    seed: u64,
) -> Result<SolveReport, String> {
    if inst.koblitz.is_none() {
        return Err("ic.shared_rank needs a Koblitz instance".into());
    }
    let (name, params) = Params::parse_spec(&m.params["factor_base"])?;
    if name != "compact-orbit-scan"
        || params.0.len() != 2
        || !params.0.contains_key("columns")
        || !params.0.contains_key("raw_x_cap")
    {
        return Err("ic.shared_rank needs exactly compact-orbit-scan:columns=N,raw_x_cap=M".into());
    }
    let rank_spec = SharedRankSpec {
        columns: usize::try_from(params.u64("columns")?).map_err(|_| "columns do not fit usize")?,
        raw_x_cap: params.u64("raw_x_cap")?,
        rank_seed: param_u64(m, "rank_seed")?,
        max_trials: param_u64(m, "rank_max_trials")?,
    };
    let target_spec = SharedTargetSpec {
        residual_seed: seed,
        max_attempts: u32::try_from(param_u64(m, "target_max_attempts")?)
            .map_err(|_| "target_max_attempts do not fit u32")?,
    };
    let calib = pinned_only_calibration("koblitz", &curve.slug);
    let mut point_bytes = b"ecbench-shared-rank-point-v1\0".to_vec();
    point_bytes.extend_from_slice(&target.x.to_le_bytes());
    point_bytes.extend_from_slice(&target.y.to_le_bytes());
    let started = Instant::now();
    let report = run_shared_rank_targets(
        inst,
        &rank_spec,
        &target_spec,
        &[target],
        sha256_hex(&point_bytes),
        &calib,
    )?;
    let wall = started.elapsed().as_nanos() as u64;
    let row = report
        .targets
        .first()
        .ok_or("shared rank returned no target")?;
    let recovered = row.recovered_log;
    let phases = vec![
        phase_of("factor_base", &report.rank.base),
        phase_of("oracle_setup", &report.rank.table),
        phase_of("relations", &report.rank.rank_search),
        phase_of("linear_algebra", &report.rank.linear_algebra),
        phase_of("verify", &report.rank.verification),
        phase_of("target_query", &row.query),
        phase_of("target_PDP", &row.pdp),
        phase_of("target_relation_check", &row.relation_check),
        phase_of("target_descent", &row.descent),
        phase_of("target_recovery_check", &row.recovery_check),
    ];
    let mut counters = BTreeMap::new();
    counters.insert("rank_trials".into(), report.rank.trials);
    counters.insert("rank_hits".into(), report.rank.hits);
    counters.insert("matrix_rank".into(), report.rank.rank as u64);
    counters.insert("target_attempts".into(), row.attempts.len() as u64);
    let mut unpriced = Vec::new();
    for phase in &phases {
        for (name, &count) in &phase.native {
            if count == 0 {
                continue;
            }
            let priced = PRICED_NATIVE
                .iter()
                .find(|(counter, _)| *counter == name)
                .is_some_and(|(_, unit)| calib.is_pinned(unit));
            if !priced {
                let label = if name.ends_with("_uncharged") {
                    name.clone()
                } else {
                    format!("{name}_uncharged")
                };
                *counters.entry(label.clone()).or_insert(0) += count;
                unpriced.push(label);
            }
        }
    }
    // The current group-addition-equivalent calibration does not account
    // for these operations, even if every named native counter is pinned.
    unpriced.extend([
        "field_arithmetic_uncharged".into(),
        "hash_and_allocation_uncharged".into(),
        "modular_combination_uncharged".into(),
    ]);
    unpriced.sort();
    unpriced.dedup();
    let builder = CompactOrbitScanBase { instance: inst };
    let group = BinaryGroup(&inst.fast);
    let ctx = InstanceCtx {
        group: &group,
        generator: inst.generator,
        target: group.identity(),
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: curve.slug.clone(),
        field_degree: Some(inst.n),
    };
    let mut scratch = GroupOps::default();
    let fb = builder.build(&ctx, &params, &mut scratch)?;
    let facts = factor_base_facts(&curve.slug, &m.params["factor_base"], &params, &fb, |p| {
        group.key(p)
    })?;
    if facts.columns != report.rank.columns as u64
        || facts.signed_points != report.rank.signed_points as u64
    {
        return Err("shared rank factor-base inventory differs from solver setup".into());
    }
    let online_names = [
        ("target_query", row.query.wall_ns),
        ("target_PDP", row.pdp.wall_ns),
        ("target_relation_check", row.relation_check.wall_ns),
        ("target_descent", row.descent.wall_ns),
        ("target_recovery_check", row.recovery_check.wall_ns),
    ];
    let online = OnlineWindow {
        wall_ns: row.online_wall_ns,
        phases_ns: online_names
            .iter()
            .map(|(name, ns)| ((*name).into(), *ns))
            .collect(),
        zero_phases: online_names
            .iter()
            .filter(|(_, ns)| *ns == 0)
            .map(|(name, _)| ((*name).into(), "phase performed no work".into()))
            .collect(),
        start_event: "Q subgroup validation after target-blind base, table and rank logs are ready"
            .into(),
        stop_event: "[d]G = Q verified by full-point scalar replay".into(),
        included_stages: online_names.iter().map(|(name, _)| (*name).into()).collect(),
        mapping: "direct Q, then run-seeded Q+[a]G residuals; exact m3 folded-table PDP, full-point witness check, column-log combination and full-point scalar replay".into(),
    };
    if online.phases_ns.values().sum::<u64>() != online.wall_ns {
        return Err("shared rank online phases do not sum to online wall interval".into());
    }
    Ok(SolveReport {
        recovered,
        exhausted: !report.verified,
        phases,
        total_gae: report.total_gae,
        automorphisms_used: 2 * inst.n,
        counters,
        unpriced,
        deterministic: true,
        nondeterminism: Vec::new(),
        solve_wall_ns: wall,
        factor_base: Some(facts),
        detail: serde_json::to_value(&report).map_err(|e| e.to_string())?,
        online: Some(online),
        online_error: None,
        solver: None,
    })
}

#[cfg(test)]
mod shared_rank_tests {
    use super::*;
    use crate::cryptanalysis::ecbench::workload::CurveSpec;

    #[test]
    fn shared_rank_adapter_keeps_the_one_target_window_and_inventory() {
        let spec = CurveSpec::Koblitz { a: 1, n: 17 };
        let inst = spec.build().unwrap();
        let facts = inst.facts(&spec);
        let Instance::Binary(binary) = &inst else {
            unreachable!()
        };
        let target = binary.fast.mul_u64(binary.generator, 113);
        let method = resolve(&MethodSpec {
            id: "ic.shared_rank".into(),
            params: [
                (
                    "factor_base".into(),
                    "compact-orbit-scan:columns=4,raw_x_cap=10000".into(),
                ),
                ("rank_seed".into(), "7".into()),
                ("rank_max_trials".into(), "10000".into()),
                ("target_max_attempts".into(), "64".into()),
            ]
            .into_iter()
            .collect(),
        })
        .unwrap();
        let report = solve(
            &method,
            &inst,
            &facts,
            &[format!("0x{:x}", target.x), format!("0x{:x}", target.y)],
            19,
        )
        .unwrap();
        assert_eq!(report.recovered, Some(113 % binary.r));
        assert_eq!(report.factor_base.as_ref().unwrap().columns, 4);
        assert!(report
            .unpriced
            .contains(&"field_arithmetic_uncharged".into()));
        let online = report
            .online
            .expect("library online clock survives ecbench");
        assert!(online.wall_ns > 0);
        assert_eq!(online.phases_ns.len(), 5);
        assert_eq!(online.phases_ns.values().sum::<u64>(), online.wall_ns);
    }
}

// ── Factor-base dumps ──────────────────────────────────────────────

pub const FB_DUMP_SCHEMA: &str = "ecbench.factor_base_dump/v1";

/// One point of a dumped factor base.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FbPoint {
    pub x: String,
    pub y: String,
    /// The relation-matrix column the point contributes to.
    pub col: u64,
    /// Its coefficient in that column (mod r).
    pub coef: u64,
}

/// A factor base with its points: what `ecbench fb` writes and the
/// database's `factor_base_points` table holds.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FactorBaseDump {
    pub schema: String,
    pub curve: CurveFacts,
    pub factor_base: FactorBaseFacts,
    /// Group operations the builder charged.
    pub build_adds: u64,
    pub build_doubles: u64,
    pub points: Vec<FbPoint>,
}

fn dump_points<E: Copy>(
    fb: &crate::cryptanalysis::ic_boundary::FactorBase<E>,
    xy: impl Fn(&E) -> (u64, u64),
) -> Vec<FbPoint> {
    fb.points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let (x, y) = xy(p);
            FbPoint {
                x: format!("0x{x:x}"),
                y: format!("0x{y:x}"),
                col: fb.col_of.get(i).copied().unwrap_or(0) as u64,
                coef: fb.coef_of.get(i).copied().unwrap_or(0),
            }
        })
        .collect()
}

/// Build the factor base `fb_spec` (a plug-in spec such as
/// `koblitz-orbit:divisor=1`) on `curve` and return it whole.  The id is
/// the one an `ic.pipeline` run with the same `factor_base` reports.
pub fn dump_factor_base(
    curve: &crate::cryptanalysis::ecbench::workload::CurveSpec,
    fb_spec: &str,
) -> Result<FactorBaseDump, String> {
    let inst = curve.build()?;
    let facts = inst.facts(curve);
    let (name, params) = Params::parse_spec(fb_spec)?;
    match &inst {
        Instance::Prime(i) => {
            let abscissa = PrimeAbscissaBase { instance: i };
            let orbit = GlvOrbitBase { instance: i };
            let base: &dyn FactorBaseBuilder<_> = match name.as_str() {
                "prime-abscissa" => &abscissa,
                "glv-orbit" => &orbit,
                other => {
                    return Err(format!(
                        "factor base `{other}` does not run on a prime curve"
                    ))
                }
            };
            let ctx = InstanceCtx {
                group: &i.curve,
                generator: i.generator_point(),
                target: i.generator_point(),
                r: i.r,
                cofactor: i.cofactor,
                group_order: i.group_order,
                name: facts.slug.clone(),
                field_degree: None,
            };
            let mut ops = GroupOps::default();
            let fb = base.build(&ctx, &params, &mut ops)?;
            let fbf = factor_base_facts(&facts.slug, fb_spec, &params, &fb, |p| i.curve.key(p))?;
            let total = {
                let mut t = fb.cost.group_ops;
                t.merge(ops);
                t
            };
            Ok(FactorBaseDump {
                schema: FB_DUMP_SCHEMA.into(),
                points: dump_points(&fb, |p| (p.x, p.y)),
                curve: facts,
                factor_base: fbf,
                build_adds: total.adds,
                build_doubles: total.doubles,
            })
        }
        Instance::Binary(i) => {
            let g = BinaryGroup(&i.fast);
            let subspace = BinarySubspaceBase { instance: i };
            let orbit = KoblitzOrbitBase { instance: i };
            let compact = CompactOrbitScanBase { instance: i };
            let base: &dyn FactorBaseBuilder<BinaryGroup> = match name.as_str() {
                "binary-subspace" => &subspace,
                "koblitz-orbit" => &orbit,
                "compact-orbit-scan" => &compact,
                other => {
                    return Err(format!(
                        "factor base `{other}` does not run on a binary curve"
                    ))
                }
            };
            let ctx = InstanceCtx {
                group: &g,
                generator: i.generator,
                target: i.generator,
                r: i.r,
                cofactor: i.cofactor,
                group_order: i.group_order,
                name: facts.slug.clone(),
                field_degree: Some(i.n),
            };
            let mut ops = GroupOps::default();
            let fb = base.build(&ctx, &params, &mut ops)?;
            let fbf = factor_base_facts(&facts.slug, fb_spec, &params, &fb, |p| g.key(p))?;
            let total = {
                let mut t = fb.cost.group_ops;
                t.merge(ops);
                t
            };
            Ok(FactorBaseDump {
                schema: FB_DUMP_SCHEMA.into(),
                points: dump_points(&fb, |p| (p.x, p.y)),
                curve: facts,
                factor_base: fbf,
                build_adds: total.adds,
                build_doubles: total.doubles,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_written_before_hashing() {
        let a = resolve(&MethodSpec {
            id: "rho.negation".into(),
            params: BTreeMap::new(),
        })
        .unwrap();
        let mut p = BTreeMap::new();
        p.insert("cap_multiple".to_string(), "64".to_string());
        let b = resolve(&MethodSpec {
            id: "rho.negation".into(),
            params: p,
        })
        .unwrap();
        assert_eq!(a.method_id, b.method_id);
        assert!(a.method_id.starts_with("ECM1h"));
    }

    #[test]
    fn solver_stats_read_the_sat_and_macaulay_figures_by_their_names() {
        use crate::cryptanalysis::ic_framework::stages::SystemShape;
        use crate::cryptanalysis::ic_framework::SolverReport;
        let shape = SystemShape {
            n_vars: 10,
            n_equations: 13,
            degrees: vec![2; 13],
            semi_regular_degree: Some(3),
        };
        let mut extra = BTreeMap::new();
        extra.insert("variables_max".to_string(), 132);
        extra.insert("clauses_max".to_string(), 2157);
        extra.insert("decisions".to_string(), 53_360);
        extra.insert("propagations".to_string(), 2_214_306);
        let sat = SolverReport {
            name: "sat-cdcl".into(),
            calls: 179,
            ops: 48_872,
            op_unit: "conflicts".into(),
            extra,
            ..Default::default()
        };
        let s = SolverStats::from_report(&sat, Some(&shape));
        assert_eq!(s.sat_conflicts, Some(48_872));
        assert_eq!(s.sat_variables, Some(132));
        assert_eq!(s.sat_clauses, Some(2157));
        assert_eq!(s.sat_decisions, Some(53_360));
        assert_eq!(
            (s.n_vars, s.n_equations, s.semi_regular_degree),
            (Some(10), Some(13), Some(3))
        );
        assert_eq!(s.macaulay_rows, None);
        // A round trip through the record's JSON keeps every field.
        let back: SolverStats = serde_json::from_str(&serde_json::to_string(&s).unwrap()).unwrap();
        assert_eq!(back, s);

        let mut extra = BTreeMap::new();
        extra.insert("matrix_rows_max".to_string(), 4096);
        extra.insert("matrix_cols_max".to_string(), 8192);
        extra.insert("max_poly_degree".to_string(), 4);
        let f4 = SolverReport {
            name: "f4-f2".into(),
            calls: 3,
            ops: 10,
            op_unit: "word XORs".into(),
            solving_degree_max: 3,
            extra,
            ..Default::default()
        };
        let s = SolverStats::from_report(&f4, None);
        assert_eq!(
            (s.macaulay_rows, s.macaulay_columns, s.macaulay_degree),
            (Some(4096), Some(8192), Some(4))
        );
        assert_eq!(s.sat_conflicts, None);
        assert_eq!(s.n_vars, None);
    }

    #[test]
    fn unknown_and_missing_parameters_are_refused() {
        let mut p = BTreeMap::new();
        p.insert("jumps".to_string(), "8".to_string());
        assert!(resolve(&MethodSpec {
            id: "bsgs.textbook".into(),
            params: p
        })
        .is_err());
        assert!(resolve(&MethodSpec {
            id: "ic.pipeline".into(),
            params: BTreeMap::new()
        })
        .is_err());
        assert!(resolve(&MethodSpec {
            id: "rho.nope".into(),
            params: BTreeMap::new()
        })
        .is_err());
    }
}
