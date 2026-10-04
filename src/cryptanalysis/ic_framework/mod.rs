//! # A pluggable index-calculus benchmarking framework.
//!
//! The repository's boundary ledger measures a fixed list of variants
//! end to end.  That is the right shape for *reporting* a result and
//! the wrong shape for *exploring* one: adding a variant means editing
//! the ladder, and comparing someone else's Gröbner engine against
//! ours means them forking it.
//!
//! This framework turns the same pipeline inside out.  Each stage is a
//! plug point, a run is a **choice at each stage**, and the report
//! carries every stage's cost in one unit so the effect of a choice is
//! visible where it lands rather than only in the total.
//!
//! ```text
//!   instance ──► factor base ──► targets ──► decomposition ──► relations ──► logarithm
//!                    │              │             │                │             │
//!            FactorBaseBuilder  Targets   DecompositionOracle  RelationSolver  verified
//!                                               │
//!                                         SystemSolver          ← F4, F5, XL, SAT, …
//! ```
//!
//! ## What a run answers
//!
//! The question the framework exists for is not "how fast is index
//! calculus" but **"what does changing *this* do to everything
//! downstream"**.  A report therefore carries, for one configuration:
//!
//! - the factor base's size, column count and build cost;
//! - the decomposition stage's trials, hit rate and cost per call, and,
//!   when the oracle is algebraic, the system's shape, the degree the
//!   solver reached and the degree a structureless system would have
//!   reached;
//! - the relation matrix's rows, rank and work;
//! - whether the recovered logarithm is right;
//! - and `S`, the whole thing divided by `√r`, so two configurations on
//!   different-sized instances are comparable.
//!
//! Change the factor base and the hit rate moves, which moves the
//! trials, which moves the matrix.  The report shows all four.
//!
//! ## The rules a run obeys
//!
//! These are `AGENTS.md`'s, and the framework enforces what it can:
//!
//! - **Every phase is inside the total.**  Base construction, table
//!   setup, failed decompositions, linear algebra and verification all
//!   land in `S`.  A number that prices one phase is a *stage
//!   diagnostic* and the report labels it as one.
//! - **Operation counts are the metric.**  Wall time rides along as a
//!   practicality note and is never the headline.
//! - **A row without a verified answer is not a result.**  The runner
//!   checks the recovered logarithm against the planted one and says so
//!   in the report; a configuration that does not recover it is
//!   reported, not hidden.
//! - **Budget exhaustion is not "no solution".**  A solver that runs
//!   out of budget is counted separately from one that refuted the
//!   system, because conflating them silently lowers the hit rate.
//!
//! ## Where to go next
//!
//! `docs/ic/FRAMEWORK.md` is the guide: the stage contracts in prose, a
//! worked example of adding a solver, how to run a sweep, and how to
//! report what comes out.

pub mod linalg;
pub mod plugins;
pub mod shared_rank;
pub mod solvers;
pub mod stages;

use std::time::Instant;

use serde::Serialize;

use crate::cryptanalysis::ic_boundary::{
    collect_and_solve_with_completion, price_phase, Calibration, CompletionRule, CountedGroup,
    FactorBase, GroupOps, PhaseCost, PipelineOutcome, RelationSolver, RestartPool, TargetSource,
};
use crate::cryptanalysis::ic_measurement as measurement;
use linalg::matrix_by_name;
use stages::{DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SystemShape, Targets};

/// One configuration: a choice at each stage, with its parameters.
#[derive(Clone, Debug, Default, Serialize)]
pub struct PipelineSpec {
    pub factor_base: String,
    pub factor_base_params: Params,
    pub oracle: String,
    pub oracle_params: Params,
    /// The polynomial solver, for oracles that use one.
    pub solver: Option<String>,
    pub solver_params: Params,
    pub targets: Targets,
    /// The matrix and stop rule: `incremental-gauss`,
    /// `incremental-gauss-full-rank`, its `-checked` variant, or
    /// `structured-gauss`.
    pub linalg: String,
    /// Cap on targets drawn before the run gives up.
    pub max_trials: u64,
    pub seed: u64,
    /// Per-solver-call wall-clock budget, in seconds; zero for none.
    pub solver_budget_seconds: u64,
}

impl PipelineSpec {
    /// A one-line name for the configuration, for table rows.
    pub fn label(&self) -> String {
        let p = |params: &Params| {
            if params.0.is_empty() {
                String::new()
            } else {
                let kv: Vec<String> = params.0.iter().map(|(k, v)| format!("{k}={v}")).collect();
                format!("[{}]", kv.join(","))
            }
        };
        format!(
            "{}{} + {}{} + {}{} + {}",
            self.factor_base,
            p(&self.factor_base_params),
            self.oracle,
            p(&self.oracle_params),
            self.solver.as_deref().unwrap_or("—"),
            p(&self.solver_params),
            self.targets.name(),
        ) + &if self.linalg.is_empty() || self.linalg == "incremental-gauss" {
            String::new()
        } else {
            format!(" + {}", self.linalg)
        }
    }
}

/// What the factor-base stage cost and produced.
#[derive(Clone, Debug, Serialize)]
pub struct FactorBaseReport {
    pub name: String,
    pub description: String,
    /// Signed points in the base: both `P` and `−P` are entries.
    pub signed_points: usize,
    /// Distinct abscissae.
    pub abscissae: usize,
    /// Unknowns in the relation matrix, which is what the fold changes.
    pub columns: usize,
    /// Signed points per column: the fold, measured rather than assumed.
    pub points_per_column: f64,
    pub cost: PhaseCost,
}

/// What the decomposition stage did.
#[derive(Clone, Debug, Serialize)]
pub struct DecompositionReport {
    pub name: String,
    pub summands: u32,
    pub description: String,
    pub targets_tried: u64,
    pub relations_found: u64,
    /// Relations per target: the number the factor base moves and every
    /// downstream stage inherits.
    pub hit_rate: f64,
    /// What `prepare` cost — a pair table, or nothing — as its own
    /// phase inside `S`, so the same base reads the same beside every
    /// oracle.
    pub setup: PhaseCost,
    pub cost: PhaseCost,
    /// The algebraic system, when the oracle built one.
    pub system: Option<SystemShape>,
    /// Totals over the solver calls, when the oracle used a solver.
    pub solver: Option<SolverReport>,
}

/// What the polynomial solver cost across a run.
#[derive(Clone, Debug, Default, Serialize)]
pub struct SolverReport {
    pub name: String,
    pub calls: u64,
    pub ops: u64,
    pub op_unit: String,
    /// This host's measured nanoseconds per `op_unit`: the conversion
    /// factor a `measured` price rests on, recorded so the count can be
    /// re-priced on another host or pinned later.
    pub ns_per_op: Option<f64>,
    pub wall_ns: u64,
    pub peak_bytes: u64,
    /// Calls that ran out of budget, counted apart from refutations.
    pub budget_exceeded: u64,
    /// Group-addition equivalents the solver work was priced at, and
    /// how: `pinned` (a repository ratio) or `measured` (this host's
    /// wall time over its addition time, which §12 of the ledger note
    /// explains is not comparable across hosts).
    pub gae: f64,
    pub priced_by: String,
    pub solving_degree_mean: Option<f64>,
    pub solving_degree_max: u32,
    /// The degree a semi-regular system of the same shape would reach.
    pub semi_regular_degree: Option<u32>,
    /// Mean solving degree over that bound: below one is structure the
    /// solver exploited.
    pub degree_over_bound: Option<f64>,
    pub extra: std::collections::BTreeMap<String, u64>,
}

/// What the relation matrix did.
#[derive(Clone, Debug, Serialize)]
pub struct LinearAlgebraReport {
    pub name: String,
    pub rows: u64,
    pub rank: u64,
    pub dependent: u64,
    pub work: u64,
    pub work_unit: String,
    pub cost: PhaseCost,
}

/// One configuration, run end to end on one instance.
#[derive(Clone, Debug, Serialize)]
pub struct RunReport {
    pub instance: String,
    pub log2_r: f64,
    pub r: u64,
    pub group_order: u64,
    pub spec: PipelineSpec,
    pub label: String,
    pub factor_base: FactorBaseReport,
    pub decomposition: DecompositionReport,
    pub linear_algebra: LinearAlgebraReport,
    pub verify: PhaseCost,
    /// The recovered logarithm, and whether it is the planted one.
    pub recovered: Option<u64>,
    pub verified: bool,
    /// Group-addition equivalents over the whole pipeline.
    pub total_gae: f64,
    /// `total_gae / sqrt(r)`: the unit every row of this repository is
    /// comparable in.
    pub s: f64,
    /// `S` of a counted Pollard rho on the same instance, when the
    /// caller supplies one; the reference every row is read against.
    pub rho_s: Option<f64>,
    pub s_over_rho: Option<f64>,
    pub wall_seconds: f64,
    /// Set when the run did not recover the logarithm within
    /// `max_trials`.  A configuration that fails is reported, not
    /// dropped: the sweep column that matters most is often which
    /// choices fail.
    pub exhausted: bool,
}

/// Certify the solved logarithm of one nonzero-coefficient point in each
/// factor-base column. A relation uses `coef_of[i] * x_col` for point `P_i`;
/// after cofactor projection its full-point identity must be
/// `[coef_of[i] * x_col]G = [h]P_i`. Both scalar multiplications are charged.
/// This is independent of the relation matrix's assertion that it has rank.
fn check_full_rank_base_columns<G: CountedGroup>(
    ctx: &InstanceCtx<G>,
    fb: &FactorBase<G::Elt>,
    matrix: &dyn RelationSolver,
    cost: &mut PhaseCost,
) -> bool {
    let started = Instant::now();
    let mut representative = vec![None; fb.columns];
    for (i, (&col, &coef)) in fb.col_of.iter().zip(&fb.coef_of).enumerate() {
        if col < fb.columns && coef % ctx.r != 0 && representative[col].is_none() {
            representative[col] = Some(i);
        }
    }
    let mut valid = true;
    for (col, index) in representative.into_iter().enumerate() {
        let Some((i, x)) = index.zip(matrix.pinned(col)) else {
            cost.count("base_log_missing_columns", 1);
            valid = false;
            continue;
        };
        let exponent = ((fb.coef_of[i] as u128 * x as u128) % ctx.r as u128) as u64;
        let expected = ctx.group.mul(&mut cost.group_ops, ctx.generator, exponent);
        let projected = ctx
            .group
            .mul(&mut cost.group_ops, fb.points[i], ctx.cofactor);
        cost.count("base_log_columns_checked", 1);
        if expected != projected {
            cost.count("base_log_verification_failures", 1);
            valid = false;
        }
    }
    cost.wall_ns += started.elapsed().as_nanos() as u64;
    valid
}

/// Run one configuration end to end.
///
/// `planted` is the logarithm the runner checks the answer against.  No
/// stage receives it: a stage that could read it could fake a relation.
#[allow(clippy::too_many_arguments)]
pub fn run_pipeline<G: CountedGroup>(
    ctx: &InstanceCtx<G>,
    spec: &PipelineSpec,
    base_builder: &dyn FactorBaseBuilder<G>,
    oracle: &mut dyn DecompositionOracle<G>,
    planted: u64,
    calib: &Calibration,
    rho_s: Option<f64>,
) -> Result<RunReport, String> {
    let started = Instant::now();

    // ── Stage 1: the factor base ───────────────────────────────────
    //
    // The phase marks below are no-ops unless the caller holds an
    // `ic_measurement::Session`; with one, the factor base and the
    // oracle's tables are target-independent set-up and the online
    // interval opens at the first target-dependent query (see
    // `collect_and_solve_with`) and closes after the verified recovery.
    measurement::mark(measurement::Phase::FactorBase);
    let mut fb_ops = GroupOps::default();
    let fb_started = Instant::now();
    let fb = base_builder.build(ctx, &spec.factor_base_params, &mut fb_ops)?;
    let mut fb_cost = fb.cost.clone();
    fb_cost.group_ops.merge(fb_ops);
    fb_cost.wall_ns = fb_cost.wall_ns.max(fb_started.elapsed().as_nanos() as u64);

    // ── Stage 2: whatever the oracle precomputes ───────────────────
    //
    // Its own phase, not folded into the base: a pair table is the
    // oracle's cost, and charging it to the base made the same base
    // read 1,830 GAE beside a table oracle and 0 beside an algebraic
    // one in the first frozen solver sweep.
    measurement::mark(measurement::Phase::Precompute);
    let mut prep_ops = GroupOps::default();
    let prep_started = Instant::now();
    oracle.prepare(ctx, &fb, &spec.oracle_params, &mut prep_ops)?;
    let mut prep_cost = PhaseCost::default();
    prep_cost.group_ops.merge(prep_ops);
    prep_cost.native = oracle.setup_native();
    prep_cost.wall_ns = prep_started.elapsed().as_nanos() as u64;

    // ── Stages 3 and 4: relations, then the matrix ─────────────────
    //
    // The loop is `ic_boundary`'s, not a second copy: the target guard,
    // the collision guards and the verification are the audited ones,
    // so a framework row is comparable with a ledger row rather than
    // merely similar to it.
    let linalg_name = if spec.linalg.is_empty() {
        "incremental-gauss"
    } else {
        spec.linalg.as_str()
    };
    let mut matrix: Box<dyn RelationSolver> = matrix_by_name(linalg_name, fb.columns + 1, ctx.r)?;
    let source = match spec.targets {
        Targets::Random => TargetSource::Random,
        Targets::Walk => TargetSource::Walk,
    };
    measurement::begin_online(measurement::Phase::TargetQuery);
    let checked_full_rank = linalg_name == "incremental-gauss-full-rank-checked";
    let completion = if linalg_name == "incremental-gauss-full-rank" || checked_full_rank {
        CompletionRule::FullRank
    } else {
        CompletionRule::TargetPinned
    };
    let mut outcome: PipelineOutcome = collect_and_solve_with_completion(
        ctx.group,
        ctx.generator,
        ctx.target,
        ctx.r,
        ctx.cofactor,
        &fb,
        spec.seed,
        spec.max_trials,
        source,
        RestartPool::Lazy,
        matrix.as_mut(),
        completion,
        |ops, counters, point| oracle.decompose(ctx, &fb, ops, counters, point),
    );
    // Kept opt-in so the previously recorded full-rank sessions retain
    // byte-identical replay and charging. A checked run cannot claim a
    // verified full-rank solution solely from the target-point check.
    if checked_full_rank && outcome.verified {
        measurement::mark(measurement::Phase::RecoveryCheck);
        outcome.verified =
            check_full_rank_base_columns(ctx, &fb, matrix.as_ref(), &mut outcome.verify);
    }
    measurement::end_online();

    // ── The report ─────────────────────────────────────────────────
    let mut fb_phase = fb_cost;
    let mut rel = outcome.relations;
    let mut la = outcome.linear_algebra;
    let mut ver = outcome.verify;
    for phase in [&mut fb_phase, &mut prep_cost, &mut rel, &mut la, &mut ver] {
        price_phase(phase, calib);
    }
    // The solver's work is engine work in its own unit.  Only a unit the
    // calibration carries a ratio for is priced by count — `word XORs`,
    // the dense Macaulay row operation §5 of the ledger note priced
    // matrix-F4 in, at the pinned `ns_per_word_xor`.  Every other unit
    // prices at this host's measured wall time over its addition time
    // and says so: the first frozen sweep showed why, since a
    // Buchberger "monomial operation" costs about 50 ns here against
    // 0.4 ns for a word XOR, and pricing it at the word ratio would have
    // flattered the engine a hundredfold.  `ns_per_op` records the
    // measured conversion so the count can be re-priced later.
    let solver_report = oracle.solver_totals().map(|t| {
        let wall_per_op = (t.ops > 0).then(|| t.wall_ns as f64 / t.ops as f64);
        // `pinned` means the ratio came from the repository's table, and
        // only `Calibration::pin` can say so: a word-XOR factor that is
        // merely present was measured on this host (every bench
        // calibration measures one), and a freshly generated curve has
        // no table entry.  A count priced at the measured ratio is still
        // `measured`, with that ratio recorded as `ns_per_op`.
        let (gae, priced_by, ns_per_op) = match (t.op_unit.as_str(), calib.ns_per_word_xor) {
            ("word XORs", Some(ns)) if calib.ns_per_add > 0.0 => (
                t.ops as f64 * ns / calib.ns_per_add,
                if calib.is_pinned("ns_per_word_xor") {
                    "pinned"
                } else {
                    "measured"
                },
                Some(ns),
            ),
            _ if calib.ns_per_add > 0.0 => {
                (t.wall_ns as f64 / calib.ns_per_add, "measured", wall_per_op)
            }
            _ => (0.0, "unpriced", wall_per_op),
        };
        rel.count("solver_calls", t.calls);
        rel.count("solver_ops", t.ops);
        rel.count("solver_budget_exceeded", t.budget_exceeded);
        rel.gae += gae;
        SolverReport {
            name: t.solver.clone(),
            calls: t.calls,
            ops: t.ops,
            op_unit: t.op_unit.clone(),
            ns_per_op,
            wall_ns: t.wall_ns,
            peak_bytes: t.peak_bytes,
            budget_exceeded: t.budget_exceeded,
            gae,
            priced_by: priced_by.into(),
            solving_degree_mean: t.solving_degree_mean(),
            solving_degree_max: t.solving_degree_max,
            semi_regular_degree: t.semi_regular_degree(),
            degree_over_bound: t.degree_over_bound(),
            extra: t.extra.clone(),
        }
    });
    let total_gae = fb_phase.gae + prep_cost.gae + rel.gae + la.gae + ver.gae;
    let sqrt_r = (ctx.r as f64).sqrt();
    let s = total_gae / sqrt_r;

    let (work, work_unit) = matrix.work();

    Ok(RunReport {
        instance: ctx.name.clone(),
        log2_r: (ctx.r as f64).log2(),
        r: ctx.r,
        group_order: ctx.group_order,
        spec: spec.clone(),
        label: spec.label(),
        factor_base: FactorBaseReport {
            name: base_builder.name().to_string(),
            description: base_builder.describe(&spec.factor_base_params),
            signed_points: fb.points.len(),
            abscissae: fb.abscissae,
            columns: fb.columns,
            points_per_column: fb.points.len() as f64 / fb.columns.max(1) as f64,
            cost: fb_phase,
        },
        decomposition: DecompositionReport {
            name: oracle.name().to_string(),
            summands: oracle.summands(),
            description: oracle.describe(&spec.oracle_params),
            targets_tried: outcome.trials,
            relations_found: outcome.relations_found,
            hit_rate: outcome.relations_found as f64 / outcome.trials.max(1) as f64,
            setup: prep_cost,
            cost: rel,
            system: oracle.last_system(),
            solver: solver_report,
        },
        linear_algebra: LinearAlgebraReport {
            name: linalg_name.into(),
            rows: outcome.relations_found,
            rank: outcome.rank,
            dependent: outcome.dependent,
            work,
            work_unit: work_unit.into(),
            cost: la,
        },
        verify: ver,
        recovered: outcome.recovered,
        verified: outcome.verified && outcome.recovered == Some(planted),
        total_gae,
        s,
        rho_s,
        s_over_rho: rho_s.map(|rho| s / rho),
        wall_seconds: started.elapsed().as_secs_f64(),
        exhausted: outcome.recovered.is_none(),
    })
}

/// The table a sweep prints.
pub fn format_markdown(rows: &[RunReport]) -> String {
    let mut out = String::new();
    out.push_str(
        "| instance | configuration | base | cols | pts/col | hit rate | trials | rows | rank | D_solve | D_sr | S | vs rho | correct |\n",
    );
    out.push_str("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|:--|\n");
    for r in rows {
        let d = r.decomposition.solver.as_ref();
        let deg = d
            .and_then(|s| s.solving_degree_mean)
            .map(|v| format!("{v:.1}"))
            .unwrap_or_else(|| "—".into());
        let bound = d
            .and_then(|s| s.semi_regular_degree)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "—".into());
        let vs_rho = r
            .s_over_rho
            .map(|v| format!("{v:.2}×"))
            .unwrap_or_else(|| "—".into());
        out.push_str(&format!(
            "| {} | {} | {} | {} | {:.1} | {:.3} | {} | {} | {} | {} | {} | {:.3} | {} | {} |\n",
            r.instance,
            r.label,
            r.factor_base.signed_points,
            r.factor_base.columns,
            r.factor_base.points_per_column,
            r.decomposition.hit_rate,
            r.decomposition.targets_tried,
            r.linear_algebra.rows,
            r.linear_algebra.rank,
            deg,
            bound,
            r.s,
            vs_rho,
            if r.verified {
                "yes"
            } else if r.exhausted {
                "gave up"
            } else {
                "NO"
            },
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::plugins::{MitmOracle, PrimeAbscissaBase, SubtractOracle};
    use super::stages::{InstanceCtx, Params, Targets};
    use super::*;
    use crate::cryptanalysis::ic_boundary::{roster_prime_instance, Calibration};

    fn ctx_and_planted(
        inst: &crate::cryptanalysis::ic_boundary::PrimeInstance,
    ) -> (
        InstanceCtx<'_, crate::cryptanalysis::ic_boundary::PrimeCurve>,
        u64,
    ) {
        let g = inst.generator_point();
        let mut ops = GroupOps::default();
        let planted = 30_011 % inst.r;
        let q =
            crate::cryptanalysis::ic_boundary::CountedGroup::mul(&inst.curve, &mut ops, g, planted);
        (
            InstanceCtx {
                group: &inst.curve,
                generator: g,
                target: q,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: None,
            },
            planted,
        )
    }

    /// **A configuration must recover the planted logarithm.**  A
    /// framework that composes stages into something that does not
    /// solve the problem is a report generator, not a benchmark.
    #[test]
    fn a_composed_pipeline_recovers_the_planted_logarithm() {
        let inst = roster_prime_instance(16).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Walk,
            max_trials: 2_000_000,
            seed: 7,
            ..Default::default()
        };
        spec.factor_base_params.set("size", "24");
        spec.oracle_params.set("negation_folded", "1");
        let mut oracle = MitmOracle::new(2);
        let report = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .expect("the configuration should run");
        assert!(report.verified, "did not recover the planted logarithm");
        assert_eq!(report.recovered, Some(planted));
        assert!(report.s > 0.0 && report.s.is_finite());
        assert!(report.factor_base.columns > 0);
        assert!(report.decomposition.relations_found >= report.linear_algebra.rank);
    }

    /// **Two oracles on the same base must agree on the answer and
    /// differ on the cost.**  That is the whole proposition of the
    /// framework: swapping one stage changes what it should change and
    /// nothing else.
    #[test]
    fn swapping_the_oracle_keeps_the_answer_and_moves_the_cost() {
        let inst = roster_prime_instance(14).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "subtract".into(),
            targets: Targets::Random,
            max_trials: 2_000_000,
            seed: 5,
            ..Default::default()
        };
        spec.factor_base_params.set("size", "16");

        let mut subtract = SubtractOracle;
        let a = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut subtract,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        spec.oracle = "mitm".into();
        let mut mitm = MitmOracle::new(2);
        let b = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut mitm,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        assert!(a.verified && b.verified, "both must recover the logarithm");
        assert_eq!(
            a.recovered, b.recovered,
            "the answer cannot depend on the oracle"
        );
        assert_eq!(
            a.factor_base.columns, b.factor_base.columns,
            "the base is the same, so its column count must be"
        );
        assert_ne!(
            a.decomposition.cost.group_ops, b.decomposition.cost.group_ops,
            "two different oracles should not cost exactly the same"
        );
    }

    /// The fold is the framework's headline lever, so it is checked:
    /// folding must cut the columns without changing the points.
    #[test]
    fn the_report_shows_what_a_column_fold_does() {
        let inst = roster_prime_instance(14).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Random,
            max_trials: 2_000_000,
            seed: 3,
            ..Default::default()
        };
        spec.factor_base_params.set("size", "16");
        let mut oracle = MitmOracle::new(2);
        let r = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();
        // A prime base folds P and −P onto one column, so the base has
        // two signed points per column.
        assert!(
            (r.factor_base.points_per_column - 2.0).abs() < 1e-9,
            "expected two signed points per column, got {}",
            r.factor_base.points_per_column
        );
    }

    /// A binary instance and a planted target on it, the way `ic bench`
    /// builds them.
    fn binary_ctx_and_planted<'a>(
        inst: &'a crate::cryptanalysis::ic_boundary::BinaryInstance,
        group: &'a crate::cryptanalysis::ic_boundary::BinaryGroup<'a>,
        planted: u64,
    ) -> InstanceCtx<'a, crate::cryptanalysis::ic_boundary::BinaryGroup<'a>> {
        let mut ops = GroupOps::default();
        let q = crate::cryptanalysis::ic_boundary::CountedGroup::mul(
            group,
            &mut ops,
            inst.generator,
            planted,
        );
        InstanceCtx {
            group,
            generator: inst.generator,
            target: q,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        }
    }

    /// **The algebraic oracle's solver is priced into `S`.**  A whole
    /// run chooses its Gröbner engine, recovers the logarithm, and the
    /// solver's cost sits inside the decomposition phase rather than
    /// beside it — the difference between a stage diagnostic and a
    /// speed, which is the reason the oracle exists.
    #[test]
    fn the_algebraic_oracle_recovers_the_logarithm_and_its_solver_is_in_s() {
        use super::plugins::{BinarySubspaceBase, DescentAlgebraicOracle};
        use super::solvers::solver_by_name;
        use crate::cryptanalysis::ic_boundary::{random_binary_instance, BinaryGroup};

        let inst = random_binary_instance(11, 3, 8).expect("a curve at n = 11");
        let group = BinaryGroup(&inst.fast);
        let planted = 1 + 30_011 % (inst.r - 1);
        let ctx = binary_ctx_and_planted(&inst, &group, planted);
        let base = BinarySubspaceBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "binary-subspace".into(),
            oracle: "descent-algebraic".into(),
            solver: Some("buchberger-f2".into()),
            targets: Targets::Walk,
            max_trials: 200_000,
            seed: 7,
            ..Default::default()
        };
        spec.factor_base_params.set("dimension", "5");
        spec.oracle_params.set("m", "2");
        let mut oracle = DescentAlgebraicOracle::new(
            2,
            &inst,
            solver_by_name("buchberger-f2").unwrap(),
            Params::default(),
            Some(std::time::Duration::from_secs(60)),
        );
        // A calibration with one nanosecond per addition, so the
        // solver's measured price reads directly as its wall time: a
        // Buchberger monomial operation has no pinned ratio (it is not a
        // word XOR, and pricing it as one would flatter the engine), so
        // the conversion is this host's measured one and the report
        // says so.
        let calib = Calibration {
            ns_per_add: 1.0,
            ns_per_word_xor: Some(0.002),
            ..Default::default()
        };
        let report = run_pipeline(&ctx, &spec, &base, &mut oracle, planted, &calib, None)
            .expect("the configuration should run");

        assert!(report.verified, "did not recover the planted logarithm");
        assert_eq!(report.recovered, Some(planted));
        let solver = report
            .decomposition
            .solver
            .as_ref()
            .expect("an algebraic oracle reports its solver");
        assert_eq!(solver.name, "buchberger-f2");
        assert!(solver.calls >= 1, "the solver was never called");
        assert_eq!(solver.priced_by, "measured");
        assert!(
            solver.ops > 0 && solver.gae > 0.0,
            "solver work was not priced"
        );
        assert!(
            (solver.gae - solver.wall_ns as f64).abs() < 1e-6,
            "a measured price is the wall time over the addition time"
        );
        let ns_per_op = solver.ns_per_op.expect("a measured conversion is recorded");
        assert!(
            (ns_per_op * solver.ops as f64 - solver.wall_ns as f64).abs() < 1.0,
            "ns_per_op times ops is the wall time"
        );
        assert!(
            report.decomposition.cost.gae >= solver.gae,
            "the solver's {} GAE must sit inside the decomposition phase's {}",
            solver.gae,
            report.decomposition.cost.gae
        );
        assert!(
            report.total_gae >= report.decomposition.cost.gae,
            "the decomposition phase must sit inside the total"
        );
        assert!(
            solver.semi_regular_degree.is_some(),
            "the degree bound is derived, not optional"
        );
        assert_eq!(report.decomposition.cost.get("solver_calls"), solver.calls);
        // Every call ended one of three ways, and each way is counted:
        // a relation, a system with no liftable solution, or no
        // solution at all.  Nothing is folded into "did not decompose".
        assert!(
            report.decomposition.relations_found + oracle.unliftable <= solver.calls,
            "more outcomes than calls"
        );
        assert_eq!(
            report.decomposition.cost.get("unliftable_systems"),
            oracle.unliftable
        );

        // Without a calibration the solver is *unpriced*, and the report
        // says so instead of quietly dropping the cost from S.
        let mut oracle = DescentAlgebraicOracle::new(
            2,
            &inst,
            solver_by_name("buchberger-f2").unwrap(),
            Params::default(),
            Some(std::time::Duration::from_secs(60)),
        );
        let bare = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();
        assert_eq!(
            bare.decomposition.solver.as_ref().unwrap().priced_by,
            "unpriced"
        );
    }

    /// An engine that reports its work as `word XORs`, the one unit with
    /// a pinnable ratio: exhaustive search with its cost relabelled.
    struct WordXorEngine(Box<dyn super::stages::SystemSolver>);

    impl super::stages::SystemSolver for WordXorEngine {
        fn name(&self) -> &str {
            "word-xor-engine"
        }
        fn describe(&self) -> String {
            "exhaustive search reporting its tests as word XORs, for the pricing test".into()
        }
        fn accepts(&self, shape: &super::stages::SystemShape) -> bool {
            self.0.accepts(shape)
        }
        fn solve(
            &self,
            system: &super::stages::BooleanSystem,
            params: &Params,
            budget: Option<std::time::Duration>,
        ) -> (super::stages::SolverVerdict, super::stages::SolverCost) {
            let (verdict, mut cost) = self.0.solve(system, params, budget);
            cost.op_unit = "word XORs".into();
            (verdict, cost)
        }
    }

    /// **`pinned` means the table, and only `Calibration::pin` can say
    /// so.**  A word-XOR factor that is merely present was measured on
    /// this host, so a count priced at it is `measured` with the ratio
    /// recorded; the same count is `pinned` only once the calibration
    /// records that the table replaced the factor.
    #[test]
    fn a_word_xor_price_is_pinned_only_when_the_calibration_says_so() {
        use super::plugins::{BinarySubspaceBase, DescentAlgebraicOracle};
        use super::solvers::solver_by_name;
        use crate::cryptanalysis::ic_boundary::{random_binary_instance, BinaryGroup};

        let inst = random_binary_instance(11, 3, 8).expect("a curve at n = 11");
        let group = BinaryGroup(&inst.fast);
        let planted = 1 + 30_011 % (inst.r - 1);
        let ctx = binary_ctx_and_planted(&inst, &group, planted);
        let base = BinarySubspaceBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "binary-subspace".into(),
            oracle: "descent-algebraic".into(),
            solver: Some("word-xor-engine".into()),
            targets: Targets::Walk,
            max_trials: 200_000,
            seed: 7,
            ..Default::default()
        };
        spec.factor_base_params.set("dimension", "5");
        let engine = || {
            DescentAlgebraicOracle::new(
                2,
                &inst,
                Box::new(WordXorEngine(solver_by_name("exhaustive").unwrap())),
                Params::default(),
                None,
            )
        };

        // Measured on this host: present, not pinned.
        let mut calib = Calibration {
            ns_per_add: 2.0,
            ns_per_word_xor: Some(0.5),
            ..Default::default()
        };
        let mut oracle = engine();
        let measured =
            run_pipeline(&ctx, &spec, &base, &mut oracle, planted, &calib, None).unwrap();
        let s = measured.decomposition.solver.as_ref().unwrap();
        assert_eq!(s.op_unit, "word XORs");
        assert_eq!(s.priced_by, "measured");
        assert!(
            (s.gae - s.ops as f64 * 0.5 / 2.0).abs() < 1e-6,
            "count times the measured ratio"
        );
        assert_eq!(s.ns_per_op, Some(0.5));

        // The table replaced the factor: pinned, same arithmetic.
        calib.pinned_units = vec!["ns_per_word_xor".into()];
        let mut oracle = engine();
        let pinned = run_pipeline(&ctx, &spec, &base, &mut oracle, planted, &calib, None).unwrap();
        let p = pinned.decomposition.solver.as_ref().unwrap();
        assert_eq!(p.priced_by, "pinned");
        assert_eq!(p.ops, s.ops, "the same run does the same work");
        assert!((p.gae - s.gae).abs() < 1e-6);
        assert!(measured.verified && pinned.verified);
    }

    /// **Two oracles that solve the same decomposition problem must agree
    /// on the logarithm.**  The algebraic oracle replaces the pair table;
    /// AGENTS.md asks that a new oracle be cross-checked against the one
    /// it replaces on a full run, and this is that check on the smallest
    /// instance that exercises it.
    #[test]
    fn the_algebraic_and_combinatorial_oracles_agree_on_the_answer() {
        use super::plugins::{BinarySubspaceBase, DescentAlgebraicOracle};
        use super::solvers::solver_by_name;
        use crate::cryptanalysis::ic_boundary::{random_binary_instance, BinaryGroup};

        let inst = random_binary_instance(11, 3, 8).expect("a curve at n = 11");
        let group = BinaryGroup(&inst.fast);
        let planted = 1 + 4_099 % (inst.r - 1);
        let ctx = binary_ctx_and_planted(&inst, &group, planted);
        let base = BinarySubspaceBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "binary-subspace".into(),
            oracle: "mitm".into(),
            targets: Targets::Walk,
            max_trials: 200_000,
            seed: 11,
            ..Default::default()
        };
        spec.factor_base_params.set("dimension", "5");

        let mut mitm = MitmOracle::new(2);
        let a = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut mitm,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        spec.oracle = "descent-algebraic".into();
        spec.solver = Some("exhaustive".into());
        let mut algebraic = DescentAlgebraicOracle::new(
            2,
            &inst,
            solver_by_name("exhaustive").unwrap(),
            Params::default(),
            None,
        );
        let b = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut algebraic,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        assert!(a.verified && b.verified, "both must recover the logarithm");
        assert_eq!(
            a.recovered, b.recovered,
            "the answer cannot depend on the oracle"
        );
        assert_eq!(a.factor_base.columns, b.factor_base.columns);
        assert!(
            a.decomposition.solver.is_none(),
            "the pair table has no solver"
        );
        assert!(
            b.decomposition.solver.is_some(),
            "the algebraic oracle has one"
        );
        assert_eq!(
            b.decomposition.cost.get("unliftable_systems"),
            algebraic.unliftable
        );

        // Target by target, the two must agree on *whether* a point
        // decomposes over the base: both are complete over it, and a
        // solver solution naming a twist abscissa lifts to nothing and
        // must never be turned into a relation.
        use super::stages::{DecompositionOracle, FactorBaseBuilder};
        use crate::cryptanalysis::ic_boundary::{CountedGroup, OracleCounters};
        let mut ops = GroupOps::default();
        let fb = base
            .build(&ctx, &spec.factor_base_params, &mut ops)
            .unwrap();
        let mut mitm = MitmOracle::new(2);
        mitm.prepare(&ctx, &fb, &Params::default(), &mut ops)
            .unwrap();
        let mut algebraic = DescentAlgebraicOracle::new(
            2,
            &inst,
            solver_by_name("exhaustive").unwrap(),
            Params::default(),
            None,
        );
        algebraic
            .prepare(&ctx, &fb, &Params::default(), &mut ops)
            .unwrap();
        let (mut ctr_a, mut ctr_b) = (OracleCounters::default(), OracleCounters::default());
        let mut hits = 0u32;
        for k in 1..=150u64 {
            let point = group.mul(&mut ops, inst.generator, k);
            let pair = mitm.decompose(&ctx, &fb, &mut ops, &mut ctr_a, point);
            let alg = algebraic.decompose(&ctx, &fb, &mut ops, &mut ctr_b, point);
            assert_eq!(
                pair.is_some(),
                alg.is_some(),
                "target [{k}]G: pair table {pair:?}, algebraic {alg:?}"
            );
            if let Some(indices) = &alg {
                let mut acc = group.identity();
                for &i in indices {
                    acc = group.add(&mut ops, acc, fb.points[i]);
                }
                assert!(
                    acc == point,
                    "an algebraic decomposition that does not sum to its target"
                );
                hits += 1;
            }
        }
        assert!(
            hits > 0,
            "no target decomposed, so nothing was cross-checked"
        );
        assert_eq!(ctr_b.unliftable_systems, algebraic.unliftable);
    }

    /// **The matrix is a lever, not a fixed cost.**  Two relation
    /// matrices on the same relations must reach the same rank and the
    /// same logarithm; what may differ is the work, and that is what the
    /// `work` column is for.
    #[test]
    fn the_structured_matrix_agrees_with_the_dense_one_end_to_end() {
        let inst = roster_prime_instance(16).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Walk,
            max_trials: 2_000_000,
            seed: 7,
            linalg: "incremental-gauss".into(),
            ..Default::default()
        };
        spec.factor_base_params.set("size", "24");
        spec.oracle_params.set("negation_folded", "1");

        let mut oracle = MitmOracle::new(2);
        let dense = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        spec.linalg = "structured-gauss".into();
        let mut oracle = MitmOracle::new(2);
        let sparse = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();

        assert!(
            dense.verified && sparse.verified,
            "both must recover the logarithm"
        );
        assert_eq!(dense.recovered, sparse.recovered);
        assert_eq!(dense.linear_algebra.name, "incremental-gauss");
        assert_eq!(sparse.linear_algebra.name, "structured-gauss");
        assert_eq!(
            dense.linear_algebra.rows, sparse.linear_algebra.rows,
            "same relations"
        );
        assert_eq!(
            dense.linear_algebra.rank, sparse.linear_algebra.rank,
            "same rank"
        );
        assert_eq!(dense.linear_algebra.work_unit, "row_ops");
        assert_eq!(sparse.linear_algebra.work_unit, "row_ops");
        assert!(dense.linear_algebra.work > 0 && sparse.linear_algebra.work > 0);
        // Same relations, same oracle: the decomposition phase is
        // identical, so any difference in S is the matrix and nothing
        // else.
        assert_eq!(
            dense.decomposition.cost.group_ops, sparse.decomposition.cost.group_ops,
            "the matrix must not change what the oracle did"
        );
    }

    #[test]
    fn full_rank_matrix_name_reaches_every_column_end_to_end() {
        let inst = roster_prime_instance(16).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Walk,
            max_trials: 2_000_000,
            seed: 7,
            linalg: "incremental-gauss-full-rank".into(),
            ..Default::default()
        };
        spec.factor_base_params.set("size", "24");
        spec.oracle_params.set("negation_folded", "1");
        let mut oracle = MitmOracle::new(2);
        let result = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut oracle,
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();
        assert!(result.verified);
        assert!(!result.exhausted);
        assert_eq!(
            result.linear_algebra.rank,
            result.factor_base.columns as u64 + 1
        );
        assert_eq!(result.linear_algebra.cost.get("full_rank_reached"), 1);
        assert!(result.linear_algebra.cost.get("first_target_pinned_rank") > 0);
        assert_eq!(result.verify.get("base_log_columns_checked"), 0);
    }

    #[test]
    fn checked_full_rank_point_checks_each_column_without_changing_the_relation_stream() {
        let inst = roster_prime_instance(16).unwrap();
        let (ctx, planted) = ctx_and_planted(&inst);
        let base = PrimeAbscissaBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "prime-abscissa".into(),
            oracle: "mitm".into(),
            targets: Targets::Walk,
            max_trials: 2_000_000,
            seed: 7,
            linalg: "incremental-gauss-full-rank".into(),
            ..Default::default()
        };
        spec.factor_base_params.set("size", "24");
        spec.oracle_params.set("negation_folded", "1");
        let control = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut MitmOracle::new(2),
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();
        spec.linalg = "incremental-gauss-full-rank-checked".into();
        let checked = run_pipeline(
            &ctx,
            &spec,
            &base,
            &mut MitmOracle::new(2),
            planted,
            &Calibration::default(),
            None,
        )
        .unwrap();
        assert!(control.verified && checked.verified);
        assert_eq!(checked.recovered, control.recovered);
        assert_eq!(
            checked.decomposition.targets_tried,
            control.decomposition.targets_tried
        );
        assert_eq!(checked.linear_algebra.rank, control.linear_algebra.rank);
        assert_eq!(
            checked.decomposition.cost.group_ops,
            control.decomposition.cost.group_ops
        );
        assert_eq!(checked.verify.get("base_log_columns_checked"), 24);
        assert_eq!(checked.verify.get("base_log_missing_columns"), 0);
        assert_eq!(checked.verify.get("base_log_verification_failures"), 0);
        assert_eq!(
            checked.verify.group_ops.scalar_mults,
            control.verify.group_ops.scalar_mults + 48
        );
    }

    #[test]
    fn checked_full_rank_rejects_a_false_column_log() {
        use crate::cryptanalysis::ic_boundary::{FactorBase, IncrementalGauss, RowStatus};

        let inst = roster_prime_instance(16).unwrap();
        let (ctx, _) = ctx_and_planted(&inst);
        let neg = ctx.group.neg(ctx.generator);
        let fb = FactorBase::from_column_map(
            "one signed column".into(),
            vec![ctx.generator, neg],
            vec![0, 0],
            vec![1, ctx.r - 1],
            1,
            |p| ctx.group.key(p),
            |p| ctx.group.key(&ctx.group.neg(*p)),
            |p| p.x,
        )
        .unwrap();
        let mut matrix = IncrementalGauss::new(2, ctx.r);
        assert_eq!(
            matrix.add_row(vec![1, 0], ctx.cofactor % ctx.r),
            RowStatus::Independent
        );
        assert_eq!(matrix.add_row(vec![0, 1], 1), RowStatus::Independent);
        let mut valid_cost = PhaseCost::default();
        assert!(check_full_rank_base_columns(
            &ctx,
            &fb,
            &matrix,
            &mut valid_cost
        ));
        assert_eq!(valid_cost.get("base_log_columns_checked"), 1);

        let mut false_matrix = IncrementalGauss::new(2, ctx.r);
        false_matrix.add_row(vec![1, 0], (ctx.cofactor + 1) % ctx.r);
        false_matrix.add_row(vec![0, 1], 1);
        let mut invalid_cost = PhaseCost::default();
        assert!(!check_full_rank_base_columns(
            &ctx,
            &fb,
            &false_matrix,
            &mut invalid_cost
        ));
        assert_eq!(invalid_cost.get("base_log_verification_failures"), 1);

        let incomplete_map = FactorBase::from_column_map(
            "declared column with no representative".into(),
            vec![ctx.generator, neg],
            vec![0, 0],
            vec![1, ctx.r - 1],
            2,
            |p| ctx.group.key(p),
            |p| ctx.group.key(&ctx.group.neg(*p)),
            |p| p.x,
        )
        .unwrap();
        let mut missing_cost = PhaseCost::default();
        assert!(!check_full_rank_base_columns(
            &ctx,
            &incomplete_map,
            &matrix,
            &mut missing_cost,
        ));
        assert_eq!(missing_cost.get("base_log_missing_columns"), 1);
    }

    /// **The torsion-symmetrised configuration is a complete index
    /// calculus.**  On `K_1 / F_{2^17}` the `u`-frame base is closed
    /// under translation by the rational 2-torsion point `T = (0, 1)`
    /// (that closure is what the `w = u² + u`, `s = Σu` rewriting
    /// assumes), the symmetrised oracle's relations lift back to the
    /// base and solve for the planted logarithm, and its solver is
    /// priced inside the decomposition phase like every other oracle's.
    #[test]
    fn the_symmetrised_koblitz_pipeline_recovers_the_logarithm_over_a_torsion_closed_base() {
        use super::plugins::{KoblitzSymmetrisedBase, SymmetrisedOracle};
        use crate::cryptanalysis::ic_boundary::{koblitz_instance, BinaryGroup, CountedGroup};
        use crate::cryptanalysis::koblitz_fast::FastPoint;

        let inst = koblitz_instance(1, 17).expect("K_1 / F_{2^17} is in the roster");
        let group = BinaryGroup(&inst.fast);
        let planted = 1 + 30_011 % (inst.r - 1);
        let ctx = binary_ctx_and_planted(&inst, &group, planted);
        let base = KoblitzSymmetrisedBase { instance: &inst };
        let mut spec = PipelineSpec {
            factor_base: "koblitz-symmetrised".into(),
            oracle: "symmetrised".into(),
            targets: Targets::Walk,
            max_trials: 200_000,
            seed: 7,
            ..Default::default()
        };
        spec.factor_base_params.set("divisor", "0;1");
        spec.oracle_params.set("divisor", "0;1");
        spec.oracle_params.set("m", "2");

        // The base is closed under P ↦ P + T, and carries V's dimension.
        let fb = base
            .build(&ctx, &spec.factor_base_params, &mut GroupOps::default())
            .expect("the u-frame base builds at divisor 0;1");
        assert_eq!(fb.dimension, Some(9));
        let t = FastPoint::affine(0, 1);
        assert!(
            fb.index_of_key(group.key(&t)).is_some(),
            "T itself is in the base"
        );
        let mut ops = GroupOps::default();
        for p in &fb.points {
            let q = group.add(&mut ops, *p, t);
            if group.is_identity(&q) {
                continue; // T + T = O
            }
            assert!(
                fb.index_of_key(group.key(&q)).is_some(),
                "P + T left the base for P = ({:#x}, {:#x})",
                p.x,
                p.y
            );
        }

        let mut oracle = SymmetrisedOracle::new(2, &inst);
        let calib = Calibration {
            ns_per_add: 200.0,
            ns_per_word_xor: Some(0.4),
            ..Default::default()
        };
        let report = run_pipeline(&ctx, &spec, &base, &mut oracle, planted, &calib, None)
            .expect("the configuration should run");
        assert!(report.verified, "did not recover the planted logarithm");
        assert_eq!(report.recovered, Some(planted));
        let solver = report
            .decomposition
            .solver
            .as_ref()
            .expect("the symmetrised oracle reports its solver");
        assert!(solver.calls >= 1 && solver.ops > 0 && solver.gae > 0.0);
        assert!(
            report.decomposition.cost.gae >= solver.gae,
            "the solver is inside the phase"
        );
        assert!(report.s > 0.0 && report.s.is_finite());
    }
}
