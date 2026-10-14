//! The PDP3 decomposition of `koblitz_index_calculus` as an `ic.pipeline`
//! oracle on a Koblitz curve: the Groebner engines (inherited F4, F6-IC,
//! matrix F4, matrix F5) and the root-finding engine
//! (`engine=rootfind`, [`PairSumTable`]) on one base.
//!
//! PR #1333 measured F6-IC, which closes branches of an inherited Boolean
//! F4 search with exact factor-base geometry, in its own worker
//! (`examples/ic_tournament_worker.rs`), on its own wall clock and on one
//! curve.  This module runs the same two functions,
//! [`groebner_decompose`] and [`groebner_decompose_f6_ic`], unmodified,
//! inside the native `ic.pipeline`, so both arms are charged in ecbench's
//! counted unit beside strong rho on the same public point.
//!
//! ## The root-finding engine
//!
//! `pdp3-koblitz:m=2|3,engine=rootfind` is the decomposition
//! `koblitz_index_calculus_dlp` runs under
//! `DecompositionStrategy::PairTable` (`RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md`):
//! the pair sums `P_i + P_j` of the base are tabulated once
//! ([`PairSumTable::build_full_within`], the summand-storing tier, so a
//! hit is read and never recovered by a scan), and a target is decomposed
//! by *finding* its partner in the table rather than by row-reducing a
//! Weil-descended Semaev system.  It takes no `degree` and no
//! `node_budget` and refuses them: there is no Groebner search to bound.
//!
//! - **Charged:** the table's `|F|(|F|+1)/2` pair additions
//!   (`|F|` of them doublings), counted from the base's size since the
//!   batched builder carries no ledger; for `m = 3` every `R − P_k` on the
//!   ledger until the first hit; and the re-addition of every hit in the
//!   group before it is returned (fail-closed, as the strategy's own
//!   `decompose` re-checks).
//! - **Counted, not charged:** one `lookup` per table probe, under the
//!   `lookups` counter the other table oracles use (priced only where the
//!   curve has a pinned ratio); the table's entry count in `oracle_setup`.
//!
//! It is exact and complete for `m ∈ {2, 3}`, so against a Groebner
//! engine with an unexhausted node budget it issues the same queries and
//! finds the same relations; the comparison is one of engines on one
//! query stream, as the F4/F6 arms are.  Its `S` is a complete count of
//! group additions, while a Groebner arm's `S` is a lower bound (the
//! solver's word XORs are counted, not charged), so the two are compared
//! on the resource vector and on wall time, not by dividing one `S` by the
//! other.
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
use num_traits::One;

use crypto_lib::binary_ecc::{BinaryCurve, F2mElement};

use crypto_lib::cryptanalysis::ecbench::methods::BinaryPlugins;
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_factor_base, lift_abscissae, BinaryGroup, BinaryInstance, ColumnFold, CountedGroup,
    FactorBase, GroupOps, OracleCounters,
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
    FrobeniusFactorBase, KoblitzCurve, PairSumTable,
};

/// The unit string of the solver's count.  It matches the symmetrised and
/// inherited-F4 adapters, so the runner treats all three alike: counted,
/// reported, and left out of `S`.
pub const WORD_XOR_UNIT: &str =
    "word XORs (elimination, specialisation and linear elimination only)";

/// Environment prefixes the Groebner code reads for tuning overrides.
const TUNING_PREFIXES: [&str; 3] = ["KIC_", "F4_", "SOLVER_"];

/// A generic binary curve `y² + xy = x³ + ax² + b` (`a ∈ {0, 1}`, `b`
/// anywhere in `F_{2^n}`) seen through the [`KoblitzCurve`] type the
/// decomposers take, with the 2-power Frobenius switched off: the same
/// shape [`KoblitzCurve::isogenous_model`] gives a floor curve
/// (`frobenius_is_endomorphism = false`, `λ = 1`, `k = 1`), but with the
/// curve's own counted order, subgroup and generator rather than a Koblitz
/// class's.  Nothing the decomposition reads depends on Frobenius: the
/// Semaev system takes only `b`, the base is folded by negation alone, and
/// F6-IC's gate needs the subspace, base membership and group addition.
fn generic_view(inst: &BinaryInstance) -> Result<KoblitzCurve, String> {
    if inst.a > 1 {
        return Err(format!(
            "the generic view takes a ∈ {{0, 1}}, this curve has a = {:#x}",
            inst.a
        ));
    }
    let n = inst.n;
    let curve = BinaryCurve {
        m: n,
        irreducible: inst.irreducible.clone(),
        a: inst.gf.to_element(inst.a),
        b: inst.gf.to_element(inst.b),
        generator: inst.fast.lower(inst.generator),
        order: BigUint::from(inst.r),
        cofactor: BigUint::from(inst.cofactor),
    };
    let trace = (1i128 << n) + 1 - i128::from(inst.group_order);
    Ok(KoblitzCurve {
        a: inst.a as u8,
        n,
        curve,
        trace: trace as i64,
        group_order: BigUint::from(inst.group_order),
        subgroup_order: BigUint::from(inst.r),
        cofactor: BigUint::from(inst.cofactor),
        lambda: BigUint::one(),
        frobenius_is_endomorphism: false,
        k: 1,
        q: 2,
        a_index: inst.a,
        b_index: 0,
        subfield_basis: vec![F2mElement::one(n)],
    })
}

/// The instance the decomposers read: the run's own when it is a Koblitz
/// curve, otherwise a copy that carries [`generic_view`].  The copy differs
/// from the run's instance only in that field; the run's group, points and
/// identity are untouched.
enum View<'i> {
    Koblitz(&'i BinaryInstance),
    Generic(Box<BinaryInstance>),
}

impl<'i> View<'i> {
    fn of(inst: &'i BinaryInstance) -> Result<Self, String> {
        if inst.koblitz.is_some() {
            return Ok(Self::Koblitz(inst));
        }
        Ok(Self::Generic(Box::new(BinaryInstance {
            name: inst.name.clone(),
            n: inst.n,
            irreducible: inst.irreducible.clone(),
            fast: inst.fast.clone(),
            gf: inst.gf.clone(),
            artin_schreier: inst.artin_schreier.clone(),
            a: inst.a,
            b: inst.b,
            group_order: inst.group_order,
            r: inst.r,
            cofactor: inst.cofactor,
            generator: inst.generator,
            koblitz: Some(generic_view(inst)?),
        })))
    }

    fn instance(&self) -> &BinaryInstance {
        match self {
            Self::Koblitz(i) => i,
            Self::Generic(i) => i,
        }
    }

    fn curve(&self) -> &KoblitzCurve {
        self.instance()
            .koblitz
            .as_ref()
            .expect("a view always carries a curve")
    }
}

fn standard_subspace(
    inst: &BinaryInstance,
    params: &Params,
) -> Result<FrobeniusFactorBase, String> {
    let kc = inst
        .koblitz
        .as_ref()
        .ok_or("this instance has no curve view")?;
    let d = params
        .get("dimension")
        .ok_or("missing parameter `dimension`; #1333 used 4 at n9 and 6 at n17")?;
    let d: u32 = d
        .parse()
        .map_err(|_| format!("dimension `{d}` is not a number"))?;
    build_standard_subspace_factor_base(kc, d)
}

/// `binary-standard-subspace:dimension=d`, also registered under its first
/// name, `koblitz-standard-subspace`.  On a Koblitz curve the two are the
/// same base; on any other binary curve it is built through
/// [`generic_view`].
pub struct KoblitzStandardSubspaceBase<'i> {
    name: &'static str,
    view: Result<View<'i>, String>,
}

impl<'i> KoblitzStandardSubspaceBase<'i> {
    pub fn new(name: &'static str, instance: &'i BinaryInstance) -> Self {
        Self {
            name,
            view: View::of(instance),
        }
    }
}

impl<'a> FactorBaseBuilder<BinaryGroup<'a>> for KoblitzStandardSubspaceBase<'_> {
    fn name(&self) -> &str {
        self.name
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
        let view = self.view.as_ref().map_err(Clone::clone)?;
        let frob = standard_subspace(view.instance(), params)?;
        let mut fb = koblitz_factor_base(
            view.instance(),
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

/// Which decomposer, and which reduction engine, the oracle calls.
///
/// Each engine runs under its own default splitting rule
/// (`SplitRule::Auto`): the inherited engines split on the highest free
/// variable, the matrix engines on the lowest.  So an engine arm is the
/// whole solver choice, as in #1333's three-way comparison, not one
/// reduction kernel on a shared tree.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pdp3Engine {
    /// `groebner_decompose` with inherited F4: the control.
    InheritedF4,
    /// `groebner_decompose_f6_ic`: inherited F4 with the geometric gate.
    F6Ic,
    /// `groebner_decompose` with Boolean matrix-F4, every node from scratch.
    MatrixF4,
    /// `groebner_decompose` with Boolean matrix-F5: matrix-F4's matrices
    /// less the rows the F5 criterion predicts reduce to zero.
    MatrixF5,
    /// The root-finding engine of `DecompositionStrategy::PairTable`: a
    /// [`PairSumTable`] of the base's pair sums, probed once per target
    /// (`m = 2`) or once per `R − P_k` (`m = 3`).  No Groebner search, so
    /// no `degree` and no `node_budget`.
    Rootfind,
}

impl Pdp3Engine {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "inherited-f4" => Ok(Self::InheritedF4),
            "f6-ic" => Ok(Self::F6Ic),
            "matrix-f4" => Ok(Self::MatrixF4),
            "matrix-f5" => Ok(Self::MatrixF5),
            "rootfind" => Ok(Self::Rootfind),
            other => Err(format!(
                "unknown engine `{other}`; try inherited-f4, f6-ic, matrix-f4, matrix-f5 or rootfind"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::InheritedF4 => "inherited-f4",
            Self::F6Ic => "f6-ic",
            Self::MatrixF4 => "matrix-f4",
            Self::MatrixF5 => "matrix-f5",
            Self::Rootfind => "rootfind",
        }
    }

    /// Whether the engine solves an algebraic system at all.
    pub fn is_groebner(self) -> bool {
        !matches!(self, Self::Rootfind)
    }

    /// The reduction engine at each node of the splitting tree, or
    /// `None` for the root-finding engine, which has no tree.
    fn solver(self, max_degree: u32) -> Option<SolverEngine> {
        match self {
            Self::InheritedF4 | Self::F6Ic => Some(SolverEngine::InheritedF4 { max_degree }),
            Self::MatrixF4 => Some(SolverEngine::MatrixF4 { max_degree }),
            Self::MatrixF5 => Some(SolverEngine::MatrixF5 { max_degree }),
            Self::Rootfind => None,
        }
    }
}

/// The group operations [`PairSumTable::build_full_within`] performs for
/// a base of `points` signed points: every `P_i + P_j` with `i ≤ j`, the
/// `i = j` terms being doublings.  Counted from the base's size because
/// the batched builder (one Montgomery inversion per row) carries no
/// ledger, the way `koblitz_factor_base` counts its constructor.
pub fn pair_table_build_ops(points: usize) -> (u64, u64) {
    let n = points as u64;
    (n * n.saturating_sub(1) / 2, n)
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

/// `pdp3-koblitz:m=3,engine=inherited-f4|f6-ic,degree=D,node_budget=N`,
/// or `pdp3-koblitz:m=2|3,engine=rootfind`.
pub struct Pdp3KoblitzOracle<'i> {
    name: &'static str,
    summands: u32,
    view: Result<View<'i>, String>,
    engine: Pdp3Engine,
    solver: SolverEngine,
    node_budget: usize,
    frob: Option<FrobeniusFactorBase>,
    index_of: HashMap<(BigUint, BigUint), usize>,
    st: Option<FieldStructure>,
    shape: Option<SystemShape>,
    totals: SolverTotals,
    /// The root-finding engine's table; `None` under a Groebner engine.
    pair: Option<PairSumTable>,
    /// Witnesses the solver returned that did not re-lift through the
    /// framework's sign resolution.  Must stay zero.
    pub unliftable: u64,
    /// Pair-table hits whose summands did not re-add to the probed point
    /// in the group.  Must stay zero; a nonzero count would mean the
    /// table and the run's base disagree on a point.
    pub table_mismatches: u64,
}

impl<'i> Pdp3KoblitzOracle<'i> {
    pub fn new(name: &'static str, summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            name,
            summands,
            view: View::of(instance),
            engine: Pdp3Engine::InheritedF4,
            solver: SolverEngine::InheritedF4 { max_degree: 3 },
            node_budget: 0,
            frob: None,
            index_of: HashMap::new(),
            st: None,
            shape: None,
            totals: SolverTotals::default(),
            pair: None,
            unliftable: 0,
            table_mismatches: 0,
        }
    }

    /// The root-finding decomposition: one table probe for `m = 2`, one
    /// per `R − P_k` for `m = 3`, every hit re-added in the group on the
    /// ledger before it is returned.
    fn decompose_rootfind<'a>(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let table = self.pair.as_ref()?;
        let g = ctx.group;
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        match self.summands {
            2 => {
                counters.lookups += 1;
                table.pairs_for(point, &mut pairs);
                for &(i, j) in &pairs {
                    let (i, j) = (i as usize, j as usize);
                    if g.add(ops, fb.points[i], fb.points[j]) == point {
                        return Some(vec![i, j]);
                    }
                    self.table_mismatches += 1;
                    counters.frobfold_mismatches += 1;
                }
                None
            }
            3 => {
                for k in 0..fb.points.len() {
                    let rest = g.add(ops, point, g.neg(fb.points[k]));
                    if g.is_identity(&rest) {
                        continue;
                    }
                    counters.lookups += 1;
                    table.pairs_for(rest, &mut pairs);
                    for &(i, j) in &pairs {
                        let (i, j) = (i as usize, j as usize);
                        if g.add(ops, fb.points[i], fb.points[j]) == rest {
                            let mut found = vec![i, j, k];
                            found.sort_unstable();
                            return Some(found);
                        }
                        self.table_mismatches += 1;
                        counters.frobfold_mismatches += 1;
                    }
                }
                None
            }
            _ => None,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for Pdp3KoblitzOracle<'_> {
    fn name(&self) -> &str {
        self.name
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        if params.get("engine") == Some("rootfind") {
            return format!(
                "root finding over the standard subspace: the pair-sum table of \
                 DecompositionStrategy::PairTable, probed once per target (m = 2) or once \
                 per R − P_k (m = {}); no Groebner search",
                self.summands
            );
        }
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
            (
                "m",
                "summands; #1333 ran 3 (Semaev S4); rootfind takes 2 or 3",
            ),
            (
                "engine",
                "inherited-f4 (the control), f6-ic, matrix-f4, matrix-f5, or rootfind (the \
                 pair-sum table of DecompositionStrategy::PairTable; no Groebner search)",
            ),
            (
                "degree",
                "highest Macaulay degree built before splitting; refused under rootfind",
            ),
            (
                "node_budget",
                "reduction calls before a solve gives up; exhaustion is a failed query; \
                 refused under rootfind",
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
        ops: &mut GroupOps,
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
        let view = self.view.as_ref().map_err(Clone::clone)?;
        let kc = view.curve();
        // Parameters that change the search have no default.
        let engine = params.get("engine").ok_or(
            "missing parameter `engine`; try inherited-f4, f6-ic, matrix-f4, matrix-f5 or rootfind",
        )?;
        self.engine = Pdp3Engine::parse(engine)?;
        let mut max_degree = 0u32;
        if self.engine.is_groebner() {
            let degree = params
                .get("degree")
                .ok_or("missing parameter `degree`; #1333 ran 3")?;
            max_degree = degree
                .parse()
                .map_err(|_| format!("degree `{degree}` is not a number"))?;
            if !(2..=4).contains(&max_degree) {
                return Err(format!("degree must be 2..=4, got {max_degree}"));
            }
            self.solver = self
                .engine
                .solver(max_degree)
                .expect("every Groebner engine has a reduction engine");
            let budget = params
                .get("node_budget")
                .ok_or("missing parameter `node_budget`; #1333 ran 8192 cold")?;
            self.node_budget = budget
                .parse()
                .map_err(|_| format!("node_budget `{budget}` is not a number"))?;
        } else {
            // The root-finding engine has no Groebner search to bound.  A
            // bound it ignored would be a parameter in the method's
            // identity that changes nothing, so it is refused.
            for key in ["degree", "node_budget"] {
                if params.get(key).is_some() {
                    return Err(format!(
                        "engine=rootfind performs no Groebner search and takes no `{key}`"
                    ));
                }
            }
        }
        let dimension = fb
            .dimension
            .ok_or("pdp3-koblitz needs the koblitz-standard-subspace base")?;
        let mut base_params = params.clone();
        if base_params.get("dimension").is_none() {
            base_params.set("dimension", dimension.to_string());
        }
        let frob = standard_subspace(view.instance(), &base_params)?;
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
            if fb.index_of_key(view.instance().fast.lift(p).pack()) != Some(i) {
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
        if !self.engine.is_groebner() {
            // The table of every `P_i + P_j`, in the summand-storing tier
            // so a hit is read, never recovered by a scan.  Its additions
            // are charged to the oracle set-up phase from the base's size.
            let table =
                PairSumTable::build_full_within(kc, &frob, PairSumTable::DEFAULT_BYTE_BUDGET)
                    .ok_or_else(|| {
                        format!(
                    "the pair-sum table of {} points does not fit the {} GiB budget or the field \
                     is too wide to pack a point into a word",
                    frob.points.len(),
                    PairSumTable::DEFAULT_BYTE_BUDGET >> 30
                )
                    })?;
            // The table's own lift must agree with the run's, point for
            // point, or a probe key would name the wrong entry.  The
            // group re-addition of every hit would catch a disagreement
            // too; this catches it before any target is spent.
            for (i, p) in frob.points.iter().enumerate() {
                if table.curve().lift(p) != fb.points[i] {
                    return Err(
                        "the pair-sum table lifts a base point differently from the run".into(),
                    );
                }
            }
            let (adds, doubles) = pair_table_build_ops(frob.points.len());
            ops.adds += adds;
            ops.doubles += doubles;
            self.pair = Some(table);
            self.frob = Some(frob);
            self.shape = None;
            self.st = None;
            self.index_of.clear();
            return Ok(());
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
        if !self.engine.is_groebner() {
            return self.decompose_rootfind(ctx, fb, ops, counters, point);
        }
        let (Some(frob), Some(st), Some(shape)) = (&self.frob, &self.st, &self.shape) else {
            return None;
        };
        let view = self.view.as_ref().ok()?;
        let kc = view.curve();
        let target = view.instance().fast.lower(point);
        let decompose = match self.engine {
            Pdp3Engine::F6Ic => groebner_decompose_f6_ic,
            _ => groebner_decompose,
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

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut counts = BTreeMap::new();
        if let Some(table) = &self.pair {
            counts.insert("pair_table_entries".into(), table.len() as u64);
        }
        counts
    }

    /// A table oracle solves no system: under `rootfind` both are `None`,
    /// as for `mitm`, and the record carries no `solver` block.
    fn last_system(&self) -> Option<SystemShape> {
        if !self.engine.is_groebner() {
            return None;
        }
        self.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        if !self.engine.is_groebner() {
            return None;
        }
        Some(self.totals.clone())
    }
}

/// The base's names: the first one it had, and the one that says it takes
/// any binary curve.  Both build the same base.
const BASE_NAMES: [&str; 2] = ["koblitz-standard-subspace", "binary-standard-subspace"];

/// The oracle's names, likewise.
const ORACLE_NAMES: [&str; 2] = ["pdp3-koblitz", "pdp3-binary"];

/// This module's plug-ins, for `methods::register_binary_plugins`.
pub fn plugins() -> BinaryPlugins {
    BinaryPlugins {
        base: |name, instance| {
            let name = BASE_NAMES.iter().find(|n| **n == name)?;
            Some(Box::new(KoblitzStandardSubspaceBase::new(name, instance))
                as Box<dyn FactorBaseBuilder<BinaryGroup<'_>> + '_>)
        },
        oracle: |name, summands, instance| {
            let name = ORACLE_NAMES.iter().find(|n| **n == name)?;
            Some(Box::new(Pdp3KoblitzOracle::new(name, summands, instance))
                as Box<dyn DecompositionOracle<BinaryGroup<'_>> + '_>)
        },
    }
}
