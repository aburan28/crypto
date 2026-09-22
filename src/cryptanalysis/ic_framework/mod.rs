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

pub mod plugins;
pub mod solvers;
pub mod stages;

use std::time::Instant;

use serde::Serialize;

use crate::cryptanalysis::ic_boundary::{
    collect_and_solve_with, price_phase, Calibration, CountedGroup, GroupOps, IncrementalGauss,
    PhaseCost, PipelineOutcome, RestartPool, TargetSource,
};
use plugins::SolverHarness;
use stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost, SystemShape, Targets,
};

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
        )
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
    pub wall_ns: u64,
    pub peak_bytes: u64,
    /// Calls that ran out of budget, counted apart from refutations.
    pub budget_exceeded: u64,
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
    harness: Option<&mut SolverHarness>,
    planted: u64,
    calib: &Calibration,
    rho_s: Option<f64>,
) -> Result<RunReport, String> {
    let started = Instant::now();

    // ── Stage 1: the factor base ───────────────────────────────────
    let mut fb_ops = GroupOps::default();
    let fb_started = Instant::now();
    let fb = base_builder.build(ctx, &spec.factor_base_params, &mut fb_ops)?;
    let mut fb_cost = fb.cost.clone();
    fb_cost.group_ops.merge(fb_ops);
    fb_cost.wall_ns = fb_cost.wall_ns.max(fb_started.elapsed().as_nanos() as u64);

    // ── Stage 2: whatever the oracle precomputes ───────────────────
    let mut prep_ops = GroupOps::default();
    oracle.prepare(ctx, &fb, &spec.oracle_params, &mut prep_ops)?;
    fb_cost.group_ops.merge(prep_ops);

    // ── Stages 3 and 4: relations, then the matrix ─────────────────
    //
    // The loop is `ic_boundary`'s, not a second copy: the target guard,
    // the collision guards and the verification are the audited ones,
    // so a framework row is comparable with a ledger row rather than
    // merely similar to it.
    let mut matrix = IncrementalGauss::new(fb.columns + 1, ctx.r);
    let source = match spec.targets {
        Targets::Random => TargetSource::Random,
        Targets::Walk => TargetSource::Walk,
    };
    let outcome: PipelineOutcome = collect_and_solve_with(
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
        &mut matrix,
        |ops, counters, point| oracle.decompose(ctx, &fb, ops, counters, point),
    );

    // ── The report ─────────────────────────────────────────────────
    let mut fb_phase = fb_cost;
    let mut rel = outcome.relations;
    let mut la = outcome.linear_algebra;
    let mut ver = outcome.verify;
    for phase in [&mut fb_phase, &mut rel, &mut la, &mut ver] {
        price_phase(phase, calib);
    }
    let total_gae = fb_phase.gae + rel.gae + la.gae + ver.gae;
    let sqrt_r = (ctx.r as f64).sqrt();
    let s = total_gae / sqrt_r;

    let solver_report = harness.map(|h| {
        let totals = h.cost();
        let shape = h.shape();
        SolverReport {
            name: h.solver.name().to_string(),
            calls: 0,
            ops: totals.map(|c| c.ops).unwrap_or(0),
            op_unit: totals.map(|c| c.op_unit.clone()).unwrap_or_default(),
            wall_ns: totals.map(|c| c.wall_ns).unwrap_or(0),
            peak_bytes: totals.map(|c| c.peak_bytes).unwrap_or(0),
            budget_exceeded: h.budget_exceeded,
            solving_degree_mean: totals.and_then(|c| c.solving_degree).map(|d| d as f64),
            solving_degree_max: totals.and_then(|c| c.solving_degree).unwrap_or(0),
            semi_regular_degree: shape.and_then(|s| s.semi_regular_degree),
            degree_over_bound: None,
            extra: totals.map(|c| c.extra.clone()).unwrap_or_default(),
        }
    });

    let (work, work_unit) = {
        use crate::cryptanalysis::ic_boundary::RelationSolver as _;
        matrix.work()
    };

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
            cost: rel,
            system: oracle.last_system(),
            solver: solver_report,
        },
        linear_algebra: LinearAlgebraReport {
            name: "incremental-gauss".into(),
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
        let q = crate::cryptanalysis::ic_boundary::CountedGroup::mul(
            &inst.curve, &mut ops, g, planted,
        );
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
            None,
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
            &ctx, &spec, &base, &mut subtract, None, planted,
            &Calibration::default(), None,
        )
        .unwrap();

        spec.oracle = "mitm".into();
        let mut mitm = MitmOracle::new(2);
        let b = run_pipeline(
            &ctx, &spec, &base, &mut mitm, None, planted,
            &Calibration::default(), None,
        )
        .unwrap();

        assert!(a.verified && b.verified, "both must recover the logarithm");
        assert_eq!(a.recovered, b.recovered, "the answer cannot depend on the oracle");
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
            &ctx, &spec, &base, &mut oracle, None, planted,
            &Calibration::default(), None,
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
}
