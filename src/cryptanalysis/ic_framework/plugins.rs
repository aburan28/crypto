//! # The factor bases, oracles and matrices that ship with the framework.
//!
//! Each of these wraps a component the repository already measures, so
//! a framework run reproduces a ledger row rather than approximating
//! one.  They are also the worked examples: a new plug-in looks like
//! these, and `docs/ic/FRAMEWORK.md` walks through writing one.
//!
//! ## Factor bases
//!
//! | name | regime | what it is |
//! |:--|:--|:--|
//! | `prime-abscissa` | prime | the smallest abscissae carrying a point, one column each |
//! | `binary-subspace` | binary, Koblitz | points whose abscissa lies in an `F_2`-subspace |
//! | `koblitz-orbit` | Koblitz | the same, each signed Frobenius orbit folded onto one column |
//! | `glv-orbit` | prime | the closure of the smallest abscissae under the curve's automorphism group, each orbit folded onto one column (`j = 0`: 6 points a column, `j = 1728`: 4, generic: 2) |
//! | `gls-line` | GLS over `F_{p²}` | the `ψ`-stable line `x ∈ u·s·F_p`, each `⟨−1, ψ⟩`-orbit folded onto one column (4 points a column) |
//!
//! The column map is the interesting part.  `koblitz-orbit` gives a
//! whole orbit of `2n` points one unknown, which cuts the relations
//! needed by `2n`; setting its `no_fold` parameter gives the control
//! that shows what the fold buys.  That pair is the one place this
//! ledger records an *advance* in the `AGENTS.md` §3 sense rather than
//! engineering, and it is a two-word change in a sweep file.
//!
//! ## Decomposition oracles
//!
//! | name | summands | how |
//! |:--|:--|:--|
//! | `subtract` | 2 | test `R − P ∈ F` for every signed base point; no table |
//! | `mitm` | 2 or 3 | a pair table probed once per target, or once per `R − P` |
//! | `mitm-frobenius` | 2 or 3 | the same on a Frobenius-folded table (Koblitz) |
//!
//! The pair table is where memory buys time, and `mitm`'s
//! `negation_folded` parameter halves the additions that build it.

use std::collections::BTreeMap;
use std::time::{Duration, Instant};

use super::stages::{
    BooleanSystem, DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost,
    SolverTotals, SolverVerdict, SystemShape, SystemSolver,
};
use crate::cryptanalysis::gls_fp2::{gls_line_base, Fp2Curve, Fp2Point, GlsInstance};
use crate::cryptanalysis::glv_invariant_base::{glv_orbit_base, AutomorphismGroup};
use crate::cryptanalysis::ic_boundary::{
    binary_subspace_factor_base, decompose_mitm, decompose_mitm_frobenius, koblitz_factor_base,
    prime_factor_base, BinaryGroup, BinaryInstance, ColumnFold, CountedGroup, FactorBase,
    FrobeniusPairTable, GroupOps, OracleCounters, PairTable, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_groebner::{f4_word_ops_thread, FieldStructure, SolverEngine};
use crate::cryptanalysis::koblitz_index_calculus::build_frobenius_factor_base_from_divisor;
use crate::cryptanalysis::koblitz_symmetrised::{
    build_symmetrised_factor_base, build_symmetrised_system, frobenius_view_of_symmetrised,
    symmetrised_groebner_decompose, SymmetrisedFactorBase,
};
use crate::cryptanalysis::pq_descent_symbolic::{descend, max_n_prime, SymbolicDescent};

// ── Factor bases ───────────────────────────────────────────────────

/// The smallest abscissae of a prime-field curve that carry a point,
/// one column each.
pub struct PrimeAbscissaBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for PrimeAbscissaBase<'_> {
    fn name(&self) -> &str {
        "prime-abscissa"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the {} smallest abscissae carrying a point, one column each",
            params.get("size").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "size",
            "abscissae in the base; the family optimum is about #E^(1/3)",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let size = params.u64("size")? as usize;
        if size == 0 {
            return Err("factor base size must be positive".into());
        }
        Ok(prime_factor_base(self.instance, size))
    }
}

/// Points of a binary curve whose abscissa lies in an `F_2`-subspace
/// of the field, one column per abscissa.
pub struct BinarySubspaceBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for BinarySubspaceBase<'_> {
    fn name(&self) -> &str {
        "binary-subspace"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "points whose abscissa lies in the {}-dimensional subspace spanned by {{1, z, …}}",
            params.get("dimension").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "dimension",
            "F_2-dimension of the abscissa subspace; about 2^dimension abscissae",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let dim = params.u64("dimension")? as u32;
        if dim == 0 || dim > 30 {
            return Err(format!("subspace dimension {dim} outside 1..=30"));
        }
        let basis: Vec<u64> = (0..dim).map(|k| 1u64 << k).collect();
        Ok(binary_subspace_factor_base(self.instance, &basis))
    }
}

/// A Koblitz base with each signed Frobenius orbit folded onto one
/// column: `2n` points, one unknown.
pub struct KoblitzOrbitBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzOrbitBase<'_> {
    fn name(&self) -> &str {
        "koblitz-orbit"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the Frobenius-invariant subspace spanned by divisor factors {}, {}",
            params.get("divisor").unwrap_or("?"),
            if params.flag("no_fold") {
                "one column per abscissa (the unfolded control)"
            } else {
                "each signed Frobenius orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "divisor",
                "comma-separated indices into the factorisation of the Frobenius characteristic \
                 polynomial; the invariant subspace they span becomes the abscissa set. This is \
                 the knob that chooses the base, so sweeping it sweeps factor bases",
            ),
            (
                "no_fold",
                "1 to give every abscissa its own column: the control that shows what the fold buys",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("this instance is not a Koblitz curve, so it has no Frobenius to fold by")?;
        let raw = params
            .get("divisor")
            .ok_or("missing parameter `divisor`; try 0 or 0,1")?;
        let idx: Result<Vec<usize>, String> = raw
            .split([',', ';'])
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| format!("divisor index `{s}` is not a number"))
            })
            .collect();
        let idx = idx?;
        if idx.is_empty() {
            return Err("parameter `divisor` selected no factors".into());
        }
        let frob = build_frobenius_factor_base_from_divisor(kc, &idx)
            .ok_or_else(|| format!("no invariant subspace for divisor {idx:?} on this curve"))?;
        let fold = if params.flag("no_fold") {
            ColumnFold::Abscissa
        } else {
            ColumnFold::SignedFrobeniusOrbit
        };
        koblitz_factor_base(
            self.instance,
            &frob,
            fold,
            format!(
                "invariant subspace, divisor {idx:?}, dimension {}",
                frob.ell
            ),
        )
        .ok_or_else(|| "the invariant subspace produced no usable factor base".into())
    }
}

/// The prime-field analogue of `koblitz-orbit`: the `size` smallest
/// abscissae, closed under the curve's automorphism group and folded
/// one orbit to a column.  `group` picks the group (`auto` takes what
/// the curve has); `no_fold` keeps the closed point set and folds by
/// negation only, which is the control that shows what the fold buys
/// on the same points.  See `glv_invariant_base`.
pub struct GlvOrbitBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for GlvOrbitBase<'_> {
    fn name(&self) -> &str {
        "glv-orbit"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the closure of the {} smallest abscissae under the {} automorphism group, {}",
            params.get("size").unwrap_or("?"),
            params.get("group").unwrap_or("auto"),
            if params.flag("no_fold") {
                "folded by negation only (the control)"
            } else {
                "each orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "size",
                "seed abscissae; the base is their closure under the group (3× as many on j = 0, 2× on j = 1728)",
            ),
            (
                "group",
                "auto (what the curve has), negation, j0 or j1728",
            ),
            (
                "no_fold",
                "1 to fold the same points by negation only: the control that shows what the fold buys",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let size = params.u64("size")? as usize;
        if size == 0 {
            return Err("factor base size must be positive".into());
        }
        let group = AutomorphismGroup::parse(params.get("group").unwrap_or("auto"))?;
        let (fb, _) = glv_orbit_base(self.instance, size, group, !params.flag("no_fold"))?;
        Ok(fb)
    }
}

/// The `ψ`-stable line of a GLS twist over `F_{p²}`, folded by
/// `⟨−1, ψ⟩`; `no_fold` folds the same points by negation only.  See
/// `gls_fp2`.
pub struct GlsLineBase<'i> {
    pub instance: &'i GlsInstance,
}

impl FactorBaseBuilder<Fp2Curve> for GlsLineBase<'_> {
    fn name(&self) -> &str {
        "gls-line"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "the ψ-stable line x ∈ u·s·F_p, {}",
            if params.flag("no_fold") {
                "folded by negation only (the control)"
            } else {
                "each ⟨−1, ψ⟩-orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "no_fold",
            "1 to fold the same points by negation only: the control that shows what the fold buys",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<Fp2Curve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<Fp2Point>, String> {
        let (fb, _) = gls_line_base(self.instance, !params.flag("no_fold"))?;
        Ok(fb)
    }
}

// ── Decomposition oracles ──────────────────────────────────────────

/// Test `R − P ∈ F` for every signed base point: no table, `|F|` group
/// operations a target, and the honest baseline every table oracle is
/// an improvement on.
pub struct SubtractOracle;

impl<G: CountedGroup> DecompositionOracle<G> for SubtractOracle {
    fn name(&self) -> &str {
        "subtract"
    }

    fn summands(&self) -> u32 {
        2
    }

    fn describe(&self, _params: &Params) -> String {
        "for each signed base point P, test whether R − P is in the base".into()
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        for (i, p) in fb.points.iter().enumerate() {
            let rest = ctx.group.add(ops, point, ctx.group.neg(*p));
            counters.lookups += 1;
            let key = ctx.group.key(&rest);
            if let Some(j) = fb.index_of_key(key) {
                return Some(vec![i, j]);
            }
        }
        None
    }
}

/// Meet in the middle over a pair table: memory for time, the lever
/// that makes a two-summand search `O(1)` a target instead of `O(|F|)`.
pub struct MitmOracle {
    m: u32,
    table: Option<PairTable>,
}

impl MitmOracle {
    pub fn new(m: u32) -> Self {
        Self { m, table: None }
    }
}

impl<G: CountedGroup> DecompositionOracle<G> for MitmOracle {
    fn name(&self) -> &str {
        "mitm"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "a table of pair sums probed once per target (m = 2) or once per R − P (m = 3), {}",
            if params.flag("negation_folded") {
                "built with the negation fold"
            } else {
                "built in full"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "negation_folded",
            "1 to build only one of each ±pair, halving the table's additions",
        )]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        let table = if params.flag("negation_folded") {
            PairTable::build_negation_folded(ctx.group, fb)
        } else {
            PairTable::build(ctx.group, fb)
        };
        // The table's construction is part of the pipeline's cost.
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        decompose_mitm(ctx.group, fb, table, self.m, ops, counters, point)
    }
}

/// Meet in the middle on a table folded by the signed Frobenius: the
/// same probe over a table `2n` times smaller.
pub struct FrobeniusMitmOracle<'i> {
    m: u32,
    instance: &'i BinaryInstance,
    table: Option<FrobeniusPairTable>,
}

impl<'i> FrobeniusMitmOracle<'i> {
    pub fn new(m: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            m,
            instance,
            table: None,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for FrobeniusMitmOracle<'_> {
    fn name(&self) -> &str {
        "mitm-frobenius"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, _params: &Params) -> String {
        "a pair table folded by the signed Frobenius, one entry per orbit".into()
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        _params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        let table = FrobeniusPairTable::build(self.instance, fb)
            .ok_or("this instance has no Koblitz structure to fold by")?;
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        decompose_mitm_frobenius(ctx.group, fb, table, self.m, ops, counters, point)
    }
}

// ── The bridge to the solver plug point ────────────────────────────

/// Holds a [`SystemSolver`] and the record of what it cost, so an
/// algebraic oracle can be written without repeating the bookkeeping.
///
/// This is what makes the solver plug point reachable from a whole
/// pipeline: swap `buchberger-f2` for `sat-cdcl` and the change
/// propagates through the relation phase to the recovered key, with
/// every stage still counted.
pub struct SolverHarness {
    pub solver: Box<dyn SystemSolver>,
    pub params: Params,
    pub budget: Option<Duration>,
    last_shape: Option<SystemShape>,
    last_cost: Option<SolverCost>,
    /// Targets whose solver hit its budget.  Counted, never silently
    /// folded into "did not decompose": an oracle that did that would
    /// quietly lower its own hit rate and make the relation phase it
    /// feeds unmeasurable.
    pub budget_exceeded: u64,
}

impl SolverHarness {
    pub fn new(solver: Box<dyn SystemSolver>, params: Params, budget: Option<Duration>) -> Self {
        Self {
            solver,
            params,
            budget,
            last_shape: None,
            last_cost: None,
            budget_exceeded: 0,
        }
    }

    /// Solve one descended system, recording what it cost.
    pub fn solve(&mut self, system: &BooleanSystem) -> SolverVerdict {
        let shape = system.shape();
        self.last_shape = Some(shape.clone());
        if !self.solver.accepts(&shape) {
            self.last_cost = None;
            self.budget_exceeded += 1;
            return SolverVerdict::BudgetExceeded;
        }
        let (verdict, cost) = self.solver.solve(system, &self.params, self.budget);
        self.last_cost = Some(cost);
        if matches!(verdict, SolverVerdict::BudgetExceeded) {
            self.budget_exceeded += 1;
        }
        verdict
    }

    pub fn shape(&self) -> Option<&SystemShape> {
        self.last_shape.as_ref()
    }

    pub fn cost(&self) -> Option<&SolverCost> {
        self.last_cost.as_ref()
    }
}

// ── The algebraic oracle ───────────────────────────────────────────

/// Descend the Semaev polynomial over the base's abscissa subspace and
/// hand the boolean system to the chosen [`SystemSolver`].
///
/// This is the oracle that makes the solver plug point reach the `S`
/// column: swap `buchberger-f2` for `sat-cdcl` and the change
/// propagates through the relation phase to the recovered logarithm,
/// with every stage still counted.  It requires a factor base whose
/// abscissae form an `F_2`-subspace — `binary-subspace` and
/// `koblitz-orbit` both record theirs in `FactorBase::subspace_basis`.
///
/// ## What it charges
///
/// The descent and the solver are engine work, not group operations,
/// and are reported through [`DecompositionOracle::solver_totals`] in
/// the solver's own unit; the runner prices them into `S`.  Lifting a
/// solution back to signed base points costs real group additions,
/// which are charged to `ops` like any other oracle's.
///
/// ## The cap
///
/// The descent is symbolic — the summation polynomial expanded in
/// `F_{2^n}[v]/(v² − v)` term by term — so it builds no table and is
/// capped only by the monomial mask: `m·n' ≤ 64` boolean variables,
/// `n' ≤ 32` at two summands and `n' ≤ 21` at three.  What limits a
/// row above that is the engine, which says so through `accepts`.
pub struct DescentAlgebraicOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    v_basis: Vec<u64>,
    harness: SolverHarness,
    totals: SolverTotals,
    /// Systems the solver decided satisfiable none of whose solutions
    /// lifted to base points summing to the target.  Expected, not a
    /// defect: a summation polynomial vanishes over the algebraic
    /// closure, so a solution may name abscissae whose points live on
    /// the quadratic twist (`y ∈ F_{2^{2n}}`, never in the base), and
    /// the pair table never sees those.  Reported as
    /// `unliftable_systems` beside `lift_failures` so the hit rate can
    /// be read against what the solver actually found.
    pub unliftable: u64,
}

impl<'i> DescentAlgebraicOracle<'i> {
    pub fn new(
        summands: u32,
        instance: &'i BinaryInstance,
        solver: Box<dyn SystemSolver>,
        solver_params: Params,
        budget: Option<Duration>,
    ) -> Self {
        let totals = SolverTotals {
            solver: solver.name().to_string(),
            ..Default::default()
        };
        Self {
            summands,
            instance,
            v_basis: Vec::new(),
            harness: SolverHarness::new(solver, solver_params, budget),
            totals,
            unliftable: 0,
        }
    }

    pub fn solver_name(&self) -> &str {
        self.harness.solver.name()
    }

    /// The descended system for a target abscissa; its `lift` maps a
    /// solution back to abscissae.  `prepare` must have run.
    fn descend(&self, x_r: u64) -> SymbolicDescent {
        descend(
            &self.instance.gf,
            self.instance.b,
            x_r,
            &self.v_basis,
            self.summands,
        )
        .expect("prepare checked the dimension")
    }
}

fn system_of(sys: &SymbolicDescent) -> BooleanSystem {
    BooleanSystem {
        equations: sys.equations.clone(),
        n_vars: sys.n_vars,
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for DescentAlgebraicOracle<'_> {
    fn name(&self) -> &str {
        "descent-algebraic"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "Weil-descend S_{} over the base's abscissa subspace and solve with {}",
            self.summands + 1,
            self.harness.solver.name()
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[("m", "summands, 2 (descends S3) or 3 (descends S4)")]
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        _params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "descent-algebraic takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let basis = fb.subspace_basis.as_ref().ok_or(
            "descent-algebraic needs a factor base whose abscissae form an F_2-subspace; \
             use binary-subspace or koblitz-orbit",
        )?;
        let cap = max_n_prime(self.summands);
        if basis.len() as u32 > cap {
            return Err(format!(
                "the descent's monomial mask holds n' ≤ {cap} at m = {}; this base has n' = {}",
                self.summands,
                basis.len()
            ));
        }
        self.v_basis = basis.clone();
        // Every system on one base has the same shape, so descend one
        // abscissa now and let the engine decline the shape here, once,
        // with a reason a sweep can report — rather than target by
        // target, which would run the whole relation phase for nothing.
        let shape = system_of(&self.descend(self.instance.generator.x)).shape();
        if !self.harness.solver.accepts(&shape) {
            return Err(format!(
                "solver `{}` declines a system of {} unknowns and {} equations; \
                 use a smaller base or another engine",
                self.harness.solver.name(),
                shape.n_vars,
                shape.n_equations
            ));
        }
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        if self.v_basis.is_empty() {
            return None;
        }
        // Descend, solve, and lift every solution until one sums to R.
        let sys = self.descend(point.x);
        let system = system_of(&sys);
        let shape = system.shape();
        let verdict = self.harness.solve(&system);
        let exceeded = matches!(verdict, SolverVerdict::BudgetExceeded);
        self.totals.absorb(&shape, self.harness.cost(), exceeded);
        let SolverVerdict::Solved(solutions) = verdict else {
            return None;
        };
        for v in solutions {
            let xs = sys.lift(v);
            if let Some(indices) =
                crate::cryptanalysis::ic_boundary::lift_abscissae(ctx.group, fb, ops, &xs, point)
            {
                return Some(indices);
            }
            counters.lift_failures += 1;
        }
        self.unliftable += 1;
        counters.unliftable_systems += 1;
        None
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.totals.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

// ── The symmetrised Koblitz base and oracle ────────────────────────

/// Parse a `divisor` parameter: indices into
/// [`crate::cryptanalysis::koblitz_index_calculus::all_factors_of_x_n_minus_1`],
/// separated by `;` (a `,` would split the plug-in spec itself).
fn divisor_indices(params: &Params) -> Result<Vec<usize>, String> {
    let raw = params
        .get("divisor")
        .ok_or("missing parameter `divisor`; try 0;1 (index 0 is x + 1)")?;
    let idx: Result<Vec<usize>, String> = raw
        .split([',', ';'])
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.trim()
                .parse::<usize>()
                .map_err(|_| format!("divisor index `{s}` is not a number"))
        })
        .collect();
    let idx = idx?;
    if idx.is_empty() {
        return Err("parameter `divisor` selected no factors".into());
    }
    Ok(idx)
}

/// `F_u = { P : u(P) ∈ V }` for the frame `u = 1/(x + 1)` and a
/// Frobenius-invariant subspace `V ∋ 1`, each signed Frobenius orbit
/// folded onto one column.
///
/// This is the base of `research/notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`
/// §8: the rational 2-torsion point `T = (0, 1)` acts on the `x`-line
/// as `x ↦ 1/x`, the frame turns that into `u ↦ u + 1`, and a subspace
/// containing `1` is closed under it for free.  `T` itself (`u = 1`) is
/// a base point.  The abscissae of `F_u` are *not* an `F_2`-subspace,
/// so the `descent-algebraic` oracle cannot run on this base; the
/// `symmetrised` oracle is the one that can, and `mitm-frobenius` runs
/// on it like on any other.  The same `divisor` handed to
/// `koblitz-orbit` gives the `x`-frame base `F_x` for the same `V`,
/// which is the matched control.
pub struct KoblitzSymmetrisedBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzSymmetrisedBase<'_> {
    fn name(&self) -> &str {
        "koblitz-symmetrised"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "F_u = {{P : 1/(x + 1) ∈ V}} for the Frobenius-invariant V ∋ 1 of divisor {}, {}",
            params.get("divisor").unwrap_or("?"),
            if params.flag("no_fold") {
                "one column per abscissa (the unfolded control)"
            } else {
                "each signed Frobenius orbit folded onto one column"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            (
                "divisor",
                "`;`-separated indices into the factorisation of xⁿ − 1 (index 0 is x + 1, which \
                 must be included so that 1 ∈ V); the same indices on koblitz-orbit give the \
                 x-frame base for the same V",
            ),
            (
                "no_fold",
                "1 to give every abscissa its own column: the control that shows what the fold buys",
            ),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("this instance is not a Koblitz curve, so it has no 1/(x + 1) frame")?;
        let idx = divisor_indices(params)?;
        let fb_u = build_symmetrised_factor_base(kc, &idx).ok_or_else(|| {
            format!(
                "no symmetrised base for divisor {idx:?}: the subspace must contain 1 \
                 (include index 0, the factor x + 1) and be a proper invariant subspace"
            )
        })?;
        let view = frobenius_view_of_symmetrised(kc, &fb_u)
            .ok_or("the u-frame base has no Frobenius-orbit view; this is a bug")?;
        let fold = if params.flag("no_fold") {
            ColumnFold::Abscissa
        } else {
            ColumnFold::SignedFrobeniusOrbit
        };
        let mut fb = koblitz_factor_base(
            self.instance,
            &view,
            fold,
            format!(
                "u-frame invariant subspace, divisor {idx:?}, dimension {}, {} points including T",
                fb_u.ell,
                fb_u.points.len()
            ),
        )
        .ok_or_else(|| "the u-frame subspace produced no usable factor base".to_string())?;
        // The base is the subspace V of the u-frame, so a frozen copy
        // carries its dimension and basis (single words, as the
        // x-frame bases do) rather than only a description.
        fb.dimension = Some(fb_u.ell as u32);
        fb.subspace_basis = Some(
            fb_u.v_basis
                .iter()
                .map(|v| u64::try_from(v.to_biguint()).unwrap_or(0))
                .collect(),
        );
        Ok(fb)
    }
}

/// Decompose over `F_u` by solving the symmetrised summation polynomial
/// in the invariants `w = u² + u`, `s = Σu` with the repository's
/// boolean matrix-F4 and splitting (`koblitz_symmetrised`).
///
/// The system has `m(ℓ − 1) + 1` unknowns of boolean degree 2 at
/// `m = 2` and 4 at `m = 3`, against `m·ℓ` unknowns of degree 2 and 6
/// for the direct `x`-frame descent on the same `V`.  A root lifts to
/// points summing to `R` or to `R + T`; both are relations over `F_u`,
/// and the second carries `T` as an extra summand.
///
/// ## What it charges
///
/// The solve is engine work, reported through `solver_totals` in the
/// same unit and by the same rule as the `inherited-f4` adapter (word
/// XORs counted, priced by the runner); lifting a root back to signed
/// base points goes through the framework's own `lift_abscissae`, so
/// its group additions are charged to `ops` exactly as the
/// `descent-algebraic` oracle's are.
///
/// ## Parameters
///
/// `m` (2 or 3), `divisor` (the base's; copied from the factor base
/// when omitted on the command line), `engine`
/// (`inherited-f4` default, `matrix-f4`, `matrix-f5`), `max_degree`
/// (Macaulay cap, default 3: the cap is absolute, so the degree-4
/// `m = 3` system still builds its own degree at the root) and
/// `node_budget` (splits before a call gives up, default 4096).
pub struct SymmetrisedOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    fb_u: Option<SymmetrisedFactorBase>,
    st: Option<FieldStructure>,
    engine: SolverEngine,
    node_budget: usize,
    shape: Option<SystemShape>,
    totals: SolverTotals,
    /// Solves whose root lifted in the `u`-frame but whose abscissae
    /// did not re-lift through the framework's sign resolution; must
    /// stay zero, and is reported so that it can be seen to.
    pub unliftable: u64,
    /// Targets with `u(R) ∈ {0, ∞}` or an oversize layout, which the
    /// system cannot be built for.
    pub unbuildable: u64,
}

impl<'i> SymmetrisedOracle<'i> {
    pub fn new(summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            summands,
            instance,
            fb_u: None,
            st: None,
            engine: SolverEngine::InheritedF4 { max_degree: 3 },
            node_budget: 4096,
            shape: None,
            totals: SolverTotals::default(),
            unliftable: 0,
            unbuildable: 0,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for SymmetrisedOracle<'_> {
    fn name(&self) -> &str {
        "symmetrised"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "solve the symmetrised S_{} in w = u² + u, s = Σu over F_u with {} (Macaulay cap {}, {} splits)",
            self.summands + 1,
            params.get("engine").unwrap_or("inherited-f4"),
            params.get("max_degree").unwrap_or("3"),
            params.get("node_budget").unwrap_or("4096"),
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands, 2 (symmetrised S3) or 3 (symmetrised S4)"),
            (
                "divisor",
                "the koblitz-symmetrised base's divisor; copied from the base when omitted",
            ),
            ("engine", "inherited-f4 (default), matrix-f4 or matrix-f5"),
            (
                "max_degree",
                "highest Macaulay degree built before splitting (default 3)",
            ),
            (
                "node_budget",
                "splits before a call gives up (default 4096)",
            ),
        ]
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<(), String> {
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "symmetrised takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("symmetrised needs a Koblitz instance")?;
        let idx = divisor_indices(params)?;
        let fb_u = build_symmetrised_factor_base(kc, &idx)
            .ok_or_else(|| format!("no symmetrised base for divisor {idx:?}"))?;
        // The base the runner built must be this F_u, point for point:
        // a mismatched divisor would make every relation a lie.
        if fb_u.points.len() != fb.points.len() {
            return Err(format!(
                "symmetrised needs the koblitz-symmetrised base with the same divisor: \
                 F_u has {} points, the base has {}",
                fb_u.points.len(),
                fb.points.len()
            ));
        }
        for p in &fb_u.points {
            let key = self.instance.fast.lift(p).pack();
            if fb.index_of_key(key).is_none() {
                return Err(
                    "symmetrised needs the koblitz-symmetrised base with the same divisor: \
                     a point of F_u is not in the base"
                        .into(),
                );
            }
        }
        let max_degree = params.u64_or("max_degree", 3)? as u32;
        self.engine = match params.get("engine").unwrap_or("inherited-f4") {
            "inherited-f4" => SolverEngine::InheritedF4 { max_degree },
            "matrix-f4" => SolverEngine::MatrixF4 { max_degree },
            "matrix-f5" => SolverEngine::MatrixF5 { max_degree },
            other => {
                return Err(format!(
                    "unknown engine `{other}`; try inherited-f4, matrix-f4 or matrix-f5"
                ))
            }
        };
        self.node_budget = params.u64_or("node_budget", 4096)? as usize;
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        // Every system on one base has the same shape; build one for the
        // generator so the report's degree columns come from a real
        // system, and so an oversize layout is refused here, once.
        let g = self.instance.fast.lower(self.instance.generator);
        let sys = build_symmetrised_system(kc, &fb_u, &g, self.summands as usize, &st).ok_or_else(
            || {
                format!(
                    "the symmetrised system for dimension {} at m = {} exceeds the 64-variable \
                     layout or the generator has u ∈ {{0, ∞}}",
                    fb_u.ell, self.summands
                )
            },
        )?;
        self.shape = Some(
            BooleanSystem {
                equations: sys.equations,
                n_vars: sys.n_vars,
            }
            .shape(),
        );
        self.totals = SolverTotals {
            solver: format!(
                "symmetrised {} (cap {}, budget {})",
                match self.engine {
                    SolverEngine::InheritedF4 { .. } => "inherited-f4",
                    SolverEngine::MatrixF4 { .. } => "matrix-f4",
                    SolverEngine::MatrixF5 { .. } => "matrix-f5",
                    _ => "engine",
                },
                max_degree,
                self.node_budget
            ),
            ..Default::default()
        };
        self.fb_u = Some(fb_u);
        self.st = Some(st);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let (Some(fb_u), Some(st), Some(shape)) = (&self.fb_u, &self.st, &self.shape) else {
            return None;
        };
        let kc = self.instance.koblitz.as_ref()?;
        let target = self.instance.fast.lower(point);
        let before = f4_word_ops_thread();
        let started = Instant::now();
        let outcome = symmetrised_groebner_decompose(
            kc,
            fb_u,
            st,
            &target,
            self.summands as usize,
            self.engine,
            self.node_budget,
        );
        let word_xors = f4_word_ops_thread().wrapping_sub(before);
        let wall_ns = started.elapsed().as_nanos() as u64;
        let Some(out) = outcome else {
            self.unbuildable += 1;
            return None;
        };
        let mut extra = BTreeMap::new();
        extra.insert("splits".to_string(), out.effort);
        extra.insert("oversize".to_string(), out.oversize as u64);
        extra.insert("used_t".to_string(), u64::from(out.used_t));
        let cost = SolverCost {
            ops: word_xors,
            // The same unit string as the `inherited-f4` adapter, so the
            // runner prices both arms by the same rule.
            op_unit: "word XORs (elimination, specialisation and linear elimination only)".into(),
            wall_ns,
            peak_bytes: 0,
            degree_reached: Some(out.built_degree),
            solving_degree: None,
            timed_out: !out.complete,
            extra,
        };
        self.totals.absorb(shape, Some(&cost), !out.complete);
        let indices = out.relation?;
        // Re-lift through the framework so the sign resolution is charged
        // exactly as it is for the x-frame oracle.  `T` (x = 0) is a base
        // point like any other, so a relation through `R + T` lifts too.
        let xs: Vec<u64> = indices
            .iter()
            .map(|&i| self.instance.fast.lift(&fb_u.points[i]).x)
            .collect();
        match crate::cryptanalysis::ic_boundary::lift_abscissae(ctx.group, fb, ops, &xs, point) {
            Some(found) => Some(found),
            None => {
                counters.lift_failures += 1;
                self.unliftable += 1;
                None
            }
        }
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}
