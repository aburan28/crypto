//! The PDP3 Groebner decomposition of `koblitz_index_calculus`, inherited
//! F4 or F6-IC, as an `ic.pipeline` oracle on a Koblitz curve.
//!
//! PR #1333 measured F6-IC, which closes branches of an inherited Boolean
//! F4 search with exact factor-base geometry, in its own worker
//! (`examples/ic_tournament_worker.rs`), on its own wall clock and on one
//! curve.  This module runs the same two functions,
//! [`groebner_decompose`] and [`groebner_decompose_f6_ic`], unmodified,
//! inside the native `ic.pipeline`, so both arms are charged in ecbench's
//! counted unit beside strong rho on the same public point.
//!
//! ## The base
//!
//! `koblitz-standard-subspace:dimension=d` is
//! [`build_standard_subspace_factor_base`]: every point whose abscissa lies
//! in the span of `1, z, …, z^{d−1}` in the curve's polynomial basis.  It is
//! the base #1333's F4/F5/F6 workers used.  The span is not
//! Frobenius-invariant, so its orbits are singletons and the column fold
//! is by negation only.
//!
//! ## What the oracle charges
//!
//! - Every affine point addition the decomposition performs is charged on
//!   the `GroupOps` ledger.  For F6-IC that includes the gate's
//!   `geometric_group_additions` (support closure, residual arithmetic and
//!   independent witness replay), which the solver counts but does not put
//!   on a ledger.  Lifting a witness to signed base points goes through the
//!   framework's `lift_abscissae`, as for every other algebraic oracle.
//! - The Boolean solver's own work is reported through `solver_totals` in
//!   word XORs, read from [`f4_word_ops_thread`] around each call, under the
//!   unit string the other inherited-F4 adapters use.  The runner leaves
//!   that unit unpriced, so `S` is a lower bound and the record says so.
//!   Reductions, splits, propagations and every F6 gate counter ride in the
//!   totals' `extra` map.
//!
//! The word-XOR counter is per thread.  `ic.pipeline` decomposes on its
//! calling thread, and the elimination runs there.
//!
//! ## Determinism
//!
//! The solver reads several tuning variables from the environment
//! (`SOLVER_SPLIT_RULE`, `KIC_*`, `F4_*`).  An override would change the
//! search without changing the method's identity, so `prepare` refuses to
//! run while any of them is set, as #1333's exclusive worker does.
//!
//! ## Why it lives in the binary
//!
//! The frozen-source replay workflows rebuild the library against pre-F6
//! snapshots of `koblitz_index_calculus.rs`
//! (`research/notes/ecc2k130/compact_frozen_source_replay_20260929`), which
//! lack the function this module calls. Other frozen evaluations pin
//! `Cargo.toml` by hash, and the tournament admits no root build script, so
//! neither a feature nor a `cfg` can gate it. The module therefore lives in
//! the `ecbench` binary, which those replays never build, and reaches
//! `ic.pipeline` through `methods::register_binary_plugins`.

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use num_bigint::BigUint;

use crypto_lib::cryptanalysis::ecbench::methods::BinaryPlugins;
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_factor_base, lift_abscissae, BinaryGroup, BinaryInstance, ColumnFold, FactorBase,
    GroupOps, OracleCounters,
};
use crypto_lib::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost, SolverTotals,
    SystemShape,
};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_word_ops_thread, FieldStructure, SolveStats, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, groebner_decompose_f6_ic, point_key,
    FrobeniusFactorBase,
};

/// The unit string of the solver's count.  It matches the symmetrised and
/// inherited-F4 adapters, so the runner treats all three alike: counted,
/// reported, and left out of `S`.
pub const WORD_XOR_UNIT: &str =
    "word XORs (elimination, specialisation and linear elimination only)";

/// Environment prefixes the Groebner code reads for tuning overrides.
const TUNING_PREFIXES: [&str; 3] = ["KIC_", "F4_", "SOLVER_"];

fn standard_subspace(
    inst: &BinaryInstance,
    params: &Params,
) -> Result<FrobeniusFactorBase, String> {
    let kc = inst
        .koblitz
        .as_ref()
        .ok_or("this instance is not a Koblitz curve")?;
    let d = params
        .get("dimension")
        .ok_or("missing parameter `dimension`; #1333 used 4 at n9 and 6 at n17")?;
    let d: u32 = d
        .parse()
        .map_err(|_| format!("dimension `{d}` is not a number"))?;
    build_standard_subspace_factor_base(kc, d)
}

/// `koblitz-standard-subspace:dimension=d`.
pub struct KoblitzStandardSubspaceBase<'i> {
    pub instance: &'i BinaryInstance,
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzStandardSubspaceBase<'_> {
    fn name(&self) -> &str {
        "koblitz-standard-subspace"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "every point with abscissa in span(1, z, …, z^{{d−1}}), d = {}, folded by negation \
             (the base of #1333's F4/F5/F6 workers)",
            params.get("dimension").unwrap_or("?")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "dimension",
            "d, the subspace dimension: 1..=20 and below the field degree",
        )]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<FastPoint>, String> {
        let frob = standard_subspace(self.instance, params)?;
        let mut fb = koblitz_factor_base(
            self.instance,
            &frob,
            ColumnFold::SignedFrobeniusOrbit,
            format!(
                "standard subspace, dimension {}, {} points",
                frob.ell,
                frob.points.len()
            ),
        )
        .ok_or("the standard subspace produced no usable factor base")?;
        fb.dimension = Some(frob.ell);
        fb.subspace_basis = Some((0..frob.ell).map(|i| 1u64 << i).collect());
        Ok(fb)
    }
}

/// Which of #1333's two decomposers the oracle calls.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pdp3Engine {
    /// `groebner_decompose` with inherited F4: the control.
    InheritedF4,
    /// `groebner_decompose_f6_ic`: inherited F4 with the geometric gate.
    F6Ic,
}

impl Pdp3Engine {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "inherited-f4" => Ok(Self::InheritedF4),
            "f6-ic" => Ok(Self::F6Ic),
            other => Err(format!(
                "unknown engine `{other}`; try inherited-f4 or f6-ic"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::InheritedF4 => "inherited-f4",
            Self::F6Ic => "f6-ic",
        }
    }
}

/// The per-call counters `extra` carries, by name.  `*_max` keys are
/// maxima over calls; the rest are totals (`SolverTotals::absorb`).
fn stats_extra(s: &SolveStats) -> BTreeMap<String, u64> {
    let mut e = BTreeMap::new();
    let mut put = |k: &str, v: u64| {
        e.insert(k.to_string(), v);
    };
    put("reductions", s.reductions as u64);
    put("splits", s.splits as u64);
    put("propagations", s.propagations as u64);
    put("infeasible_branches", s.infeasible_branches as u64);
    put("eliminated", s.eliminated as u64);
    put("oversize", s.oversize as u64);
    put("exhausted", u64::from(s.exhausted));
    put("unsupported", u64::from(s.unsupported));
    put("max_degree_built_max", u64::from(s.max_degree_built));
    put("geometric_refutations", s.geometric_refutations as u64);
    put("geometric_witnesses", s.geometric_witnesses as u64);
    put("geometric_support_checks", s.geometric_support_checks);
    put("geometric_residual_lookups", s.geometric_residual_lookups);
    put(
        "geometric_fast_residual_lookups",
        s.geometric_fast_residual_lookups,
    );
    put("geometric_batch_groups", s.geometric_batch_groups);
    put("geometric_group_additions", s.geometric_group_additions);
    put("geometric_fallbacks", s.geometric_fallbacks as u64);
    e
}

/// `pdp3-koblitz:m=3,engine=inherited-f4|f6-ic,degree=D,node_budget=N`.
pub struct Pdp3KoblitzOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    engine: Pdp3Engine,
    solver: SolverEngine,
    node_budget: usize,
    frob: Option<FrobeniusFactorBase>,
    index_of: HashMap<(BigUint, BigUint), usize>,
    st: Option<FieldStructure>,
    shape: Option<SystemShape>,
    totals: SolverTotals,
    /// Witnesses the solver returned that did not re-lift through the
    /// framework's sign resolution.  Must stay zero.
    pub unliftable: u64,
}

impl<'i> Pdp3KoblitzOracle<'i> {
    pub fn new(summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            summands,
            instance,
            engine: Pdp3Engine::InheritedF4,
            solver: SolverEngine::InheritedF4 { max_degree: 3 },
            node_budget: 0,
            frob: None,
            index_of: HashMap::new(),
            st: None,
            shape: None,
            totals: SolverTotals::default(),
            unliftable: 0,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for Pdp3KoblitzOracle<'_> {
    fn name(&self) -> &str {
        "pdp3-koblitz"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "Weil-descend S_{} over the standard subspace and solve with #1333's {} \
             (Macaulay cap {}, {} reduction calls)",
            self.summands + 1,
            params.get("engine").unwrap_or("?"),
            params.get("degree").unwrap_or("?"),
            params.get("node_budget").unwrap_or("?"),
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands; #1333 ran 3 (Semaev S4)"),
            ("engine", "inherited-f4 (the control) or f6-ic"),
            ("degree", "highest Macaulay degree built before splitting"),
            (
                "node_budget",
                "reduction calls before a solve gives up; exhaustion is a failed query",
            ),
            (
                "dimension",
                "the base's; copied from koblitz-standard-subspace when omitted",
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
        if let Some((k, _)) =
            std::env::vars().find(|(k, _)| TUNING_PREFIXES.iter().any(|p| k.starts_with(p)))
        {
            return Err(format!(
                "pdp3-koblitz refuses to run with the solver override `{k}` set: it would \
                 change the search without changing the method's identity"
            ));
        }
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "pdp3-koblitz takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let kc = self
            .instance
            .koblitz
            .as_ref()
            .ok_or("pdp3-koblitz needs a Koblitz instance")?;
        // Parameters that change the search have no default.
        let engine = params
            .get("engine")
            .ok_or("missing parameter `engine`; try inherited-f4 or f6-ic")?;
        self.engine = Pdp3Engine::parse(engine)?;
        let degree = params
            .get("degree")
            .ok_or("missing parameter `degree`; #1333 ran 3")?;
        let max_degree: u32 = degree
            .parse()
            .map_err(|_| format!("degree `{degree}` is not a number"))?;
        if !(2..=4).contains(&max_degree) {
            return Err(format!("degree must be 2..=4, got {max_degree}"));
        }
        self.solver = SolverEngine::InheritedF4 { max_degree };
        let budget = params
            .get("node_budget")
            .ok_or("missing parameter `node_budget`; #1333 ran 8192 cold")?;
        self.node_budget = budget
            .parse()
            .map_err(|_| format!("node_budget `{budget}` is not a number"))?;
        let dimension = fb
            .dimension
            .ok_or("pdp3-koblitz needs the koblitz-standard-subspace base")?;
        let mut base_params = params.clone();
        if base_params.get("dimension").is_none() {
            base_params.set("dimension", dimension.to_string());
        }
        let frob = standard_subspace(self.instance, &base_params)?;
        // The base the runner built must be this one, point for point and
        // in the same order: the decomposer returns indices into `frob`,
        // and the relation is written against the runner's base.
        if frob.points.len() != fb.points.len() {
            return Err(format!(
                "pdp3-koblitz needs koblitz-standard-subspace with the same dimension: its \
                 base has {} points, the run's has {}",
                frob.points.len(),
                fb.points.len()
            ));
        }
        for (i, p) in frob.points.iter().enumerate() {
            if fb.index_of_key(self.instance.fast.lift(p).pack()) != Some(i) {
                return Err(
                    "pdp3-koblitz needs koblitz-standard-subspace with the same dimension: \
                     the bases differ"
                        .into(),
                );
            }
        }
        if (self.summands as u64) * u64::from(frob.ell) > 64 {
            return Err(format!(
                "the descent has m·d = {}·{} > 64 Boolean unknowns",
                self.summands, frob.ell
            ));
        }
        self.index_of = frob
            .points
            .iter()
            .enumerate()
            .map(|(i, p)| (point_key(p), i))
            .collect();
        self.st = Some(FieldStructure::new(kc.n, &kc.curve.irreducible));
        // Every system on one base has one shape: m·d unknowns and n
        // equations.  The degree columns stay empty; the solver reports
        // the degree it built in `extra`.
        self.shape = Some(SystemShape {
            n_vars: (self.summands * frob.ell) as usize,
            n_equations: kc.n as usize,
            degrees: Vec::new(),
            semi_regular_degree: None,
        });
        self.totals = SolverTotals {
            solver: format!(
                "pdp3-koblitz {} (cap {}, budget {})",
                self.engine.name(),
                max_degree,
                self.node_budget
            ),
            ..Default::default()
        };
        self.frob = Some(frob);
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
        let (Some(frob), Some(st), Some(shape)) = (&self.frob, &self.st, &self.shape) else {
            return None;
        };
        let kc = self.instance.koblitz.as_ref()?;
        let target = self.instance.fast.lower(point);
        let decompose = match self.engine {
            Pdp3Engine::InheritedF4 => groebner_decompose,
            Pdp3Engine::F6Ic => groebner_decompose_f6_ic,
        };
        let before = f4_word_ops_thread();
        let started = Instant::now();
        let (indices, stats) = decompose(
            kc,
            frob,
            &self.index_of,
            st,
            &target,
            self.summands as usize,
            self.solver,
            self.node_budget,
        );
        let word_xors = f4_word_ops_thread().wrapping_sub(before);
        let wall_ns = started.elapsed().as_nanos() as u64;
        // The gate's point additions are group work: charge them.
        ops.adds += stats.geometric_group_additions;
        let cost = SolverCost {
            ops: word_xors,
            op_unit: WORD_XOR_UNIT.into(),
            wall_ns,
            peak_bytes: 0,
            degree_reached: Some(stats.max_degree_built),
            solving_degree: None,
            timed_out: stats.exhausted,
            extra: stats_extra(&stats),
        };
        self.totals.absorb(shape, Some(&cost), stats.exhausted);
        let indices = indices?;
        let xs: Vec<u64> = indices.iter().map(|&i| fb.points[i].x).collect();
        match lift_abscissae(ctx.group, fb, ops, &xs, point) {
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

/// This module's plug-ins, for `methods::register_binary_plugins`.
pub fn plugins() -> BinaryPlugins {
    BinaryPlugins {
        base: |name, instance| {
            (name == "koblitz-standard-subspace").then(|| {
                Box::new(KoblitzStandardSubspaceBase { instance })
                    as Box<dyn FactorBaseBuilder<BinaryGroup<'_>> + '_>
            })
        },
        oracle: |name, summands, instance| {
            (name == "pdp3-koblitz").then(|| {
                Box::new(Pdp3KoblitzOracle::new(summands, instance))
                    as Box<dyn DecompositionOracle<BinaryGroup<'_>> + '_>
            })
        },
    }
}
