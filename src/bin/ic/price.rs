//! `ic price`: the workflow's pipeline, every phase priced (ledger §20).
//!
//! `ic workflow` runs select → collect → logs → solve with resumable
//! state on disk, and its stage timers mix phases: the collect stage
//! builds the pair table on its first unit, the logs stage collects
//! further units, verifies relations and runs the linear algebra.  That
//! is right for a resumable run and wrong for a price.  This runs the
//! same calls in the same order, in memory and on one thread, and
//! charges each to its own exclusive clock:
//!
//! 1. setup: the curve and the options;
//! 2. select: the factor base and its projected columns, as the workflow
//!    charges them;
//! 3. build: the pair table, at the tier the workflow would choose;
//! 4. collect: the planned units and every extension unit, with the
//!    aim the workflow computes;
//! 5. logs_setup, verify, la: the log solver's setup, each relation
//!    batch's verification, each solve attempt;
//! 6. descent_setup, descent: the individual-log solver, then each
//!    target's descent (its recovery check included);
//! 7. verify_final: `[d]G = Q` and the known answer, per target.
//!
//! Target construction is outside every clock: the method receives `Q`.
//! Nothing is written to disk.
//!
//! **The unit** is one batched affine addition (`add_many` over 1,024
//! subgroup points) on the same curve, measured immediately before and
//! after every repetition; a repetition is converted at the mean of its
//! two.  Every repetition rebuilds everything from nothing, and the
//! native counts of every repetition must be identical.
//!
//! **The reference** is batch rho (`signed_frobenius_rho_batch`) on the
//! same targets, its counted group operations priced at a canonical step
//! (one batched addition plus `SignedFrobeniusClasses::canon`), measured
//! here against the same unit; Bailey et al.'s step is priced beside it
//! as a model.  With `--cold-rho-seed`, each target is also walked alone.
//!
//! The counts are what a run on another host re-prices; the prices are
//! this host's.  The control that the counts are the workflow's own is
//! run outside: `ic workflow` on the same parameter file must report the
//! same per-unit and per-target counts.
//!
//! **`--single-target`** prices the question AGENTS.md makes primary: one
//! previously unseen target, the index calculus against rho on that same
//! point, in one process under one resource envelope.
//!
//! - **Reusable set-up** is every phase above through `descent_setup`,
//!   on the exclusive clock as before, reported apart from the target.
//! - **The index calculus's online interval** starts at the first
//!   target-dependent operation, the descent's query, and stops when the
//!   descent returns a logarithm it has checked as `[d]G = Q`.  It is
//!   split by the measurement session's exclusive phase clocks
//!   (`ic_measurement`) into `target_query`, `target_pdp`,
//!   `target_relation_check`, `target_descent` and `recovery_check`,
//!   which sum to it exactly.
//! - **Rho** is [`ParallelRho`] on the same point: its jump table, in `G`
//!   only, is its reusable set-up; its online interval runs from lifting
//!   the target to the logarithm it has checked, split into `rho_solve`
//!   (walks and collisions) and `recovery_check`.
//! - **The replay**, `[d]G = Q` in the general arithmetic, is made for
//!   both arms after their intervals close and timed on its own.  Each
//!   arm's answer is also written as a certificate and its SHA-256, for a
//!   replay outside this process.
//!
//! Each repetition rebuilds both arms from nothing, the index calculus
//! first.  The parameter file must name exactly one target.
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_instance, rho_cap, signed_frobenius_rho_batch, BinaryGroup, BinaryInstance,
    CountedGroup, GroupOps, ParallelRho, ParallelRhoResult, RhoBatchResult, RhoClasses,
    SignedFrobeniusClasses,
};
use crypto_lib::cryptanalysis::ic_measurement::{self as measurement, Phase, Session};
use crypto_lib::cryptanalysis::koblitz_fast::{
    BatchScratch, FastPoint, FrobeniusCanon, FrobeniusPowers,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    projected_signed_orbit_count, CollectedRelation, ColumnCoverage, DecompositionStrategy,
    FactorBaseLogSolver, IndividualLogSolver, KoblitzCurve, RelationCollector, RelationWorkUnit,
};
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};

use super::experiment;
use super::workflow::{self, FactorBaseSource};

/// Additions a unit measurement times.
const UNIT_ADDS: usize = 2_000_000;
/// Canonicalisations a step-price measurement times.
const CANONS: usize = 2_000_000;
/// A repetition shorter than this is repeated `--repeats-fast` times.
const FAST_NS: u128 = 50_000_000;

/// Every clock, in the order the pipeline first reaches it.  The
/// `*_setup` clocks and `select_projection` hold the constructions the
/// workflow repeats around its work (orbit maps, coverage, solvers), so
/// that the work itself — `select`, `build`, `collect`, `la`, `descent`
/// — can be read apart from them; together they are the whole run.
const PHASES: [&str; 13] = [
    "setup",
    "select",
    "select_projection",
    "build",
    "collect_setup",
    "collect",
    "logs_setup",
    "verify",
    "la",
    "descent_setup",
    "descent",
    "verify_final",
    "other",
];

#[derive(Args, Clone, Debug)]
pub struct PriceArgs {
    /// The workflow parameter file to price (a `spec` factor base).
    #[arg(long)]
    pub params: PathBuf,
    /// Repetitions of the whole index-calculus pipeline.
    #[arg(long, default_value_t = 3)]
    pub repeats: usize,
    /// Repetitions instead, when the first takes under 50 ms.
    #[arg(long, default_value_t = 15)]
    pub repeats_fast: usize,
    /// Seed of batch rho on the same targets; omit to skip the reference.
    #[arg(long)]
    pub rho_seed: Option<u64>,
    /// Also walk each target alone, with seed `cold_rho_seed + i`.
    #[arg(long)]
    pub cold_rho_seed: Option<u64>,
    /// Rounds of the step-price measurement.
    #[arg(long, default_value_t = 5)]
    pub price_rounds: usize,
    /// Allow more than one Rayon thread.  The phases are then not
    /// single-thread prices, and the report says so.
    #[arg(long)]
    pub allow_threads: bool,
    /// Price one target online against rho on the same point (see the
    /// module notes).  Needs `--rho-seed` and a file with one target.
    #[arg(long)]
    pub single_target: bool,
    /// Lanes of the single-target rho; its default follows the floor.
    #[arg(long, requires = "single_target")]
    pub rho_lanes: Option<usize>,
    /// Distinguished-point bits of the single-target rho; its default
    /// follows the floor and the lanes.
    #[arg(long, requires = "single_target")]
    pub rho_dp_bits: Option<u32>,
}

/// Exclusive phase clocks: every nanosecond between two laps belongs to
/// the phase named at the second.
struct Clock {
    last: Instant,
    ns: BTreeMap<&'static str, u128>,
}

impl Clock {
    fn start() -> Self {
        Self {
            last: Instant::now(),
            ns: BTreeMap::new(),
        }
    }
    fn lap(&mut self, phase: &'static str) {
        let now = Instant::now();
        *self.ns.entry(phase).or_default() += now.duration_since(self.last).as_nanos();
        self.last = now;
    }
}

/// The unit's fixture: 1,024 subgroup points and a point to add to them.
struct UnitBench {
    g: FastPoint,
    batch: Vec<FastPoint>,
    scratch: BatchScratch,
    out: Vec<FastPoint>,
}

impl UnitBench {
    fn new(inst: &BinaryInstance) -> Self {
        let bg = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let batch: Vec<FastPoint> = (0..1024u64)
            .map(|i| {
                bg.mul(
                    &mut ops,
                    inst.generator,
                    1 + i.wrapping_mul(0x9E37_79B9_7F4A_7C15) % (inst.r - 1),
                )
            })
            .collect();
        let g = bg.mul(&mut ops, inst.generator, 1 + 0x5DEE_CE66 % (inst.r - 1));
        Self {
            g,
            batch,
            scratch: BatchScratch::default(),
            out: Vec::with_capacity(1024),
        }
    }
    /// Nanoseconds per batched addition.
    fn measure(&mut self, inst: &BinaryInstance) -> f64 {
        let reps = UNIT_ADDS / self.batch.len();
        let t = Instant::now();
        for _ in 0..reps {
            self.out.clear();
            inst.fast
                .add_many(self.g, &self.batch, &mut self.out, &mut self.scratch);
        }
        let ns = t.elapsed().as_nanos() as f64 / (reps * self.batch.len()) as f64;
        std::hint::black_box(self.out.last().map(|p| p.x));
        ns
    }
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    if s.is_empty() {
        return f64::NAN;
    }
    let m = s.len() / 2;
    if s.len() % 2 == 1 {
        s[m]
    } else {
        (s[m - 1] + s[m]) / 2.0
    }
}

/// One target, as the method receives it, with the answer it must give
/// (`None` for a public point, whose logarithm nobody constructed).
struct Target {
    q: BinaryPoint,
    expected: Option<BigUint>,
    record: workflow::TargetRecord,
}

/// What a pass does once the descent is ready.
#[derive(Clone, Copy)]
enum PassTargets<'t> {
    /// Every target in turn, each on the exclusive clock.
    Batch(&'t [Target]),
    /// One target's online interval on the measurement session.  A
    /// strict session fails closed on a phase boundary off this thread;
    /// only one can be open in a process, so tests use the plain one.
    Online { target: &'t Target, strict: bool },
}

/// The index calculus's online interval for one target.
struct OnlineIc {
    /// The session's online interval and its exclusive phases, in
    /// nanoseconds; a phase the descent never entered is `None`.
    wall_ns: u64,
    phases_ns: BTreeMap<&'static str, Option<u64>>,
    /// The same interval on a plain clock around the session's markers.
    outer_ns: u128,
    /// `[d]G = Q` in the general arithmetic, after the interval closed.
    replay_ns: u128,
    trials: usize,
    recovered: Option<BigUint>,
    /// A logarithm came back, the replay holds, and it is the known one
    /// where there is one.
    verified: bool,
}

/// The index calculus's online phases, as the measurement session names
/// them, in the order the claim schema lists them.
const IC_ONLINE_PHASES: [&str; 5] = [
    "target_query",
    "target_pdp",
    "target_relation_check",
    "target_descent",
    "recovery_check",
];

/// The claim schema's name for each of them (`docs/ic/boundary_targets.json`,
/// `vs_rho`): its stage and its field in milliseconds.
const IC_CLAIM_PHASES: [(&str, &str); 5] = [
    ("target_query", "T_target_query_ms"),
    ("target_PDP", "T_target_PDP_ms"),
    ("target_relation_check", "T_target_relation_check_ms"),
    ("target_descent", "T_target_descent_ms"),
    ("target_recovery_check", "T_target_recovery_check_ms"),
];

/// Rho's online phases.
const RHO_ONLINE_PHASES: [&str; 2] = ["rho_solve", "recovery_check"];

/// One pass of the pipeline: its phase clocks and its native counts.
struct Pass {
    ns: BTreeMap<&'static str, u128>,
    counts: Value,
    recovered: Vec<Option<String>>,
    all_verified: bool,
    /// In a single-target pass: the online interval, and the factor
    /// base's signed-orbit representatives for the identity records.
    online: Option<OnlineIc>,
    factor_base_orbits: Option<Vec<Value>>,
}

fn one_pass(p: &workflow::WorkflowParams, targets: &[Target]) -> Result<Pass, String> {
    one_pass_with(p, PassTargets::Batch(targets))
}

fn one_pass_with(p: &workflow::WorkflowParams, targets: PassTargets<'_>) -> Result<Pass, String> {
    let mut clock = Clock::start();
    let c = experiment::curve(
        p.curve.degree,
        p.curve.curve_a,
        p.curve.subfield,
        p.curve.curve_b,
    )?;
    let ic = experiment::with_linear_algebra(
        experiment::ic_options_with_descent(
            p.solver,
            p.summands,
            p.descent_summands,
            p.collection_window,
            p.max_trials,
            p.seed,
        ),
        p.linear_algebra.mode,
        p.linear_algebra.sparse,
    );
    let FactorBaseSource::Spec { spec } = &p.factor_base else {
        return Err("ic price takes a spec factor base; run the search separately".into());
    };
    clock.lap("setup");

    // Select, and the projection the workflow charges to it.
    let (fb, selection_cost) = experiment::materialize_with_selection_cost(&c, spec)?;
    clock.lap("select");
    let columns = projected_signed_orbit_count(&c, &fb);
    clock.lap("select_projection");
    if let Some(window) = p.collection_window {
        if window as usize >= fb.points.len() {
            return Err(format!(
                "collection_window must be smaller than the {}-point factor base",
                fb.points.len()
            ));
        }
    }

    // Build, at the tier the workflow's own probe budget selects.
    let pair = if ic.strategy == DecompositionStrategy::PairTable {
        Some(workflow::build_pair_table(&c, &fb, p).ok_or("field too wide for the pair table")?)
    } else {
        None
    };
    clock.lap("build");

    // Collect the planned units, aimed as the workflow aims them.
    let work = |unit: usize| RelationWorkUnit {
        seed: p.seed,
        start: unit as u64 * p.collection.unit_trials,
        count: p.collection.unit_trials,
    };
    let mut units: BTreeMap<usize, (Vec<CollectedRelation>, u64, u64)> = BTreeMap::new();
    {
        let collector = RelationCollector::with_pair_table(&c, &fb, &ic, pair.as_ref())
            .ok_or("factor base cannot decompose with this summand count")?;
        let mut coverage = if p.collection_aim {
            Some(ColumnCoverage::new(&c, &fb).ok_or("factor base has no projected columns")?)
        } else {
            None
        };
        clock.lap("collect_setup");
        for u in 0..p.collection.units {
            let aim = coverage
                .as_ref()
                .map(|cov| cov.missing_points())
                .filter(|pts| !pts.is_empty());
            let (relations, report) = collector.collect_aimed(work(u), aim.as_deref());
            if let Some(cov) = coverage.as_mut() {
                cov.add(&relations);
            }
            units.insert(
                u,
                (relations, report.trials as u64, report.summands_scanned),
            );
        }
    }
    clock.lap("collect");

    // Logs: verify, solve, and extend until determined.
    let merged: Vec<CollectedRelation> = units
        .values()
        .flat_map(|(rels, _, _)| rels.iter().cloned())
        .collect();
    let mut solver =
        FactorBaseLogSolver::new(&c, &fb, &ic).ok_or("factor base has no projected columns")?;
    clock.lap("logs_setup");
    solver.push(&merged);
    clock.lap("verify");
    let mut outcome = solver.try_solve();
    clock.lap("la");
    let mut extended = 0usize;
    if outcome.is_none() {
        let collector = RelationCollector::with_pair_table(&c, &fb, &ic, pair.as_ref())
            .ok_or("factor base cannot decompose with this summand count")?;
        let mut extend_coverage = if p.collection_aim {
            let mut cov =
                ColumnCoverage::new(&c, &fb).ok_or("factor base has no projected columns")?;
            cov.add(&merged);
            Some(cov)
        } else {
            None
        };
        clock.lap("collect_setup");
        while outcome.is_none() {
            let next = units.keys().max().map_or(0, |m| m + 1);
            if next >= p.collection.max_units {
                break;
            }
            let aim = extend_coverage
                .as_ref()
                .map(|cov| cov.missing_points())
                .filter(|pts| !pts.is_empty());
            let (relations, report) = collector.collect_aimed(work(next), aim.as_deref());
            if let Some(cov) = extend_coverage.as_mut() {
                cov.add(&relations);
            }
            clock.lap("collect");
            solver.push(&relations);
            clock.lap("verify");
            units.insert(
                next,
                (relations, report.trials as u64, report.summands_scanned),
            );
            extended += 1;
            outcome = solver.try_solve();
            clock.lap("la");
        }
    }
    let final_report = solver.report();
    let table = match outcome {
        Some((table, report)) if report.verified => table,
        _ => {
            return Err(format!(
                "relations from {} units did not determine every one of {} columns",
                units.len(),
                final_report.columns
            ))
        }
    };
    clock.lap("other");

    // Descent, then the final check the workflow makes.
    let descent = IndividualLogSolver::new(&c, &fb, &table, &ic, pair.as_ref())
        .ok_or("the log table does not match the base")?;
    clock.lap("descent_setup");
    let mut trials = Vec::new();
    let mut recovered = Vec::new();
    let mut all_verified = true;
    let mut online = None;
    let mut factor_base_orbits = None;
    match targets {
        PassTargets::Batch(targets) => {
            for t in targets {
                let report = descent.solve_report(&t.q);
                clock.lap("descent");
                let ok = report
                    .log
                    .as_ref()
                    .is_some_and(|d| c.mul(c.generator(), d) == t.q)
                    && t.expected
                        .as_ref()
                        .is_none_or(|known| report.log.as_ref() == Some(known));
                clock.lap("verify_final");
                all_verified &= ok;
                trials.push(report.trials);
                recovered.push(report.log.as_ref().map(ToString::to_string));
            }
        }
        PassTargets::Online { target, strict } => {
            let session = if strict {
                Session::begin_strict()
            } else {
                Session::begin()
            }?;
            let outer = Instant::now();
            measurement::begin_online(Phase::TargetQuery);
            let report = descent.solve_report(&target.q);
            measurement::end_online();
            let outer_ns = outer.elapsed().as_nanos();
            let snapshot = session.finish()?;
            clock.lap("online");
            let replay = Instant::now();
            let replayed = report
                .log
                .as_ref()
                .is_some_and(|d| d < &c.subgroup_order && c.mul(c.generator(), d) == target.q);
            let replay_ns = replay.elapsed().as_nanos();
            clock.lap("replay");
            let known = target
                .expected
                .as_ref()
                .is_none_or(|k| report.log.as_ref() == Some(k));
            let wall_ns = snapshot
                .online_wall_ns
                .ok_or("the online interval did not close")?;
            all_verified &= replayed && known;
            trials.push(report.trials);
            recovered.push(report.log.as_ref().map(ToString::to_string));
            online = Some(OnlineIc {
                wall_ns,
                phases_ns: IC_ONLINE_PHASES
                    .iter()
                    .map(|ph| (*ph, snapshot.online_phases_ns.get(ph).copied().flatten()))
                    .collect(),
                outer_ns,
                replay_ns,
                trials: report.trials,
                recovered: report.log.clone(),
                verified: replayed && known,
            });
            factor_base_orbits = Some(
                fb.signed_orbits
                    .iter()
                    .map(|orbit| point_pair(&fb.points[orbit[0]]))
                    .collect(),
            );
        }
    }

    let per_unit: Vec<Value> = units
        .iter()
        .map(|(u, (rels, tr, sc))| json!([u, rels.len(), tr, sc]))
        .collect();
    let counts = json!({
        "select": {
            "points": fb.points.len(),
            "signed_orbits": fb.signed_orbits.len(),
            "columns": columns,
            "selection_cost": selection_cost,
        },
        "build": {
            "tier": pair.as_ref().map(|t| t.tier()),
            "stored_pairs": pair.as_ref().map(|t| t.len()),
        },
        "collect": {
            "units": units.len(),
            "units_planned": p.collection.units,
            "units_extended": extended,
            "trials": units.values().map(|u| u.1).sum::<u64>(),
            "summands_scanned": units.values().map(|u| u.2).sum::<u64>(),
            "relations_collected": units.values().map(|u| u.0.len()).sum::<usize>(),
            "per_unit_index_relations_trials_summands": per_unit,
        },
        "logs": {
            "columns": final_report.columns,
            "relations_accepted": final_report.relations,
            "rejected": final_report.rejected_relations,
            "duplicates": final_report.duplicate_relations,
            "solve_attempts": final_report.solve_attempts,
            "linear_algebra": experiment::linear_algebra_json(&final_report)["sparse"].clone(),
        },
        "descent": {
            "summands": p.descent_summands.unwrap_or(p.summands),
            "trials_per_target": trials,
            "trials_total": trials.iter().sum::<usize>(),
        },
    });
    Ok(Pass {
        ns: clock.ns,
        counts,
        recovered,
        all_verified,
        online,
        factor_base_orbits,
    })
}

/// `[x, y]` as decimal strings, the encoding the identity records use.
fn point_pair(p: &BinaryPoint) -> Value {
    match p {
        BinaryPoint::Infinity => Value::Null,
        BinaryPoint::Affine { x, y } => {
            json!([x.to_biguint().to_string(), y.to_biguint().to_string()])
        }
    }
}

/// Lowercase hex SHA-256 of a JSON value's canonical form: serde_json's
/// compact output, whose object keys come out sorted.
fn canonical_sha256(value: &Value) -> String {
    let text = serde_json::to_string(value).expect("a JSON value serialises");
    hex::encode(sha256(text.as_bytes()))
}

/// The curve a target and a certificate are bound to.
fn curve_binding(c: &KoblitzCurve) -> Value {
    json!({
        "degree": c.n,
        "curve_a": c.a,
        "irreducible_low_terms": c.curve.irreducible.low_terms,
        "subgroup_order": c.subgroup_order.to_string(),
        "generator": point_pair(c.generator()),
    })
}

/// A recovered logarithm as a statement anyone can replay:
/// `[scalar]generator = target` on the bound curve.
fn certificate(c: &KoblitzCurve, arm: &str, method: &str, q: &BinaryPoint, d: &BigUint) -> Value {
    json!({
        "schema": "ic-price-single-target-replay-v1",
        "arm": arm,
        "method": method,
        "curve": curve_binding(c),
        "target": point_pair(q),
        "scalar": d.to_string(),
        "statement": "[scalar]generator = target, in the prime-order subgroup",
    })
}

/// The curve and target as the identity records want them
/// (`research/ic_candidate_tournament_20260915/oracle.py`'s fixture).
fn fixture(c: &KoblitzCurve, target: &Target) -> Value {
    json!({
        "degree": c.n,
        "curve_a": c.a,
        "irreducible": {
            "degree": c.curve.irreducible.degree,
            "low_terms": c.curve.irreducible.low_terms,
        },
        "subgroup_order": c.subgroup_order.to_string(),
        "group_order": c.group_order.to_string(),
        "cofactor": c.cofactor.to_string(),
        "lambda": c.lambda.to_string(),
        "generator": point_pair(c.generator()),
        "targets": [point_pair(&target.q)],
        "target_seeds": [target.record.public_hash_seed],
        "target_scalar_constructed": target.record.target_scalar_constructed,
    })
}

/// Rho's online interval on one target.
struct OnlineRho {
    wall_ns: u64,
    phases_ns: BTreeMap<&'static str, Option<u64>>,
    outer_ns: u128,
    replay_ns: u128,
    result: ParallelRhoResult,
    verified: bool,
}

/// Rho on `target` with `rho`'s reusable set-up already built: the
/// online interval from lifting the target to the logarithm the walk has
/// checked, then the replay in the general arithmetic.
fn rho_online(
    c: &KoblitzCurve,
    inst: &BinaryInstance,
    rho: &ParallelRho<'_>,
    target: &Target,
    strict: bool,
) -> Result<OnlineRho, String> {
    let session = if strict {
        Session::begin_strict()
    } else {
        Session::begin()
    }?;
    let outer = Instant::now();
    measurement::begin_online(Phase::RhoSolve);
    let q = inst.fast.lift(&target.q);
    let result = rho.solve(q, rho_cap(inst.r, 64.0));
    measurement::end_online();
    let outer_ns = outer.elapsed().as_nanos();
    let snapshot = session.finish()?;
    let replay = Instant::now();
    let replayed = result
        .recovered
        .is_some_and(|d| c.mul(c.generator(), &BigUint::from(d)) == target.q);
    let replay_ns = replay.elapsed().as_nanos();
    let known = target
        .expected
        .as_ref()
        .is_none_or(|k| result.recovered.map(BigUint::from).as_ref() == Some(k));
    Ok(OnlineRho {
        wall_ns: snapshot
            .online_wall_ns
            .ok_or("the online interval did not close")?,
        phases_ns: RHO_ONLINE_PHASES
            .iter()
            .map(|ph| (*ph, snapshot.online_phases_ns.get(ph).copied().flatten()))
            .collect(),
        outer_ns,
        replay_ns,
        verified: replayed && known,
        result,
    })
}

/// One repetition of a single-target price: both arms from nothing.
struct SingleRep {
    pass: Pass,
    rho_setup_ns: u128,
    rho: OnlineRho,
}

/// How the single-target rho is shaped.
#[derive(Clone, Copy)]
struct RhoShape {
    seed: u64,
    lanes: usize,
    dp_bits: u32,
}

fn single_rep(
    p: &workflow::WorkflowParams,
    c: &KoblitzCurve,
    target: &Target,
    shape: RhoShape,
    strict: bool,
) -> Result<SingleRep, String> {
    let pass = one_pass_with(p, PassTargets::Online { target, strict })?;
    let started = Instant::now();
    let inst = koblitz_instance(p.curve.curve_a, p.curve.degree)
        .ok_or("no single-word Koblitz instance for this curve")?;
    let rho = ParallelRho::new(&inst, shape.seed, shape.lanes, Some(shape.dp_bits))
        .ok_or("no signed-Frobenius classes")?;
    let rho_setup_ns = started.elapsed().as_nanos();
    let rho = rho_online(c, &inst, &rho, target, strict)?;
    Ok(SingleRep {
        pass,
        rho_setup_ns,
        rho,
    })
}

/// The reusable set-up's exclusive clocks in a single-target pass, and
/// the two after it.
const SINGLE_SETUP_PHASES: [&str; 11] = [
    "setup",
    "select",
    "select_projection",
    "build",
    "collect_setup",
    "collect",
    "logs_setup",
    "verify",
    "la",
    "other",
    "descent_setup",
];

/// The index of the median of `v` (the lower one for an even count).
fn median_index(v: &[f64]) -> usize {
    let mut order: Vec<usize> = (0..v.len()).collect();
    order.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    order[(v.len() - 1) / 2]
}

/// The processor list this process may run on, where the OS says.
fn cpus_allowed_list() -> Option<String> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status
        .lines()
        .find_map(|l| l.strip_prefix("Cpus_allowed_list:"))
        .map(|v| v.trim().to_string())
}

#[allow(clippy::too_many_arguments)]
fn run_single(
    args: &PriceArgs,
    p: &workflow::WorkflowParams,
    c: &KoblitzCurve,
    inst: &BinaryInstance,
    target: &Target,
    threads: usize,
    strict: bool,
    say: &dyn Fn(String),
) -> Result<Value, String> {
    let seed = args
        .rho_seed
        .ok_or("--single-target needs --rho-seed: rho runs on the same point")?;
    let lanes = args
        .rho_lanes
        .unwrap_or_else(|| ParallelRho::default_lanes(inst));
    let dp_bits = args
        .rho_dp_bits
        .unwrap_or_else(|| ParallelRho::default_dp_bits(inst, lanes));
    let shape = RhoShape {
        seed,
        lanes,
        dp_bits,
    };
    let r = inst.r;
    let sqrt_r = (r as f64).sqrt();
    let mut bench = UnitBench::new(inst);
    bench.measure(inst);

    let mut reps: Vec<Value> = Vec::new();
    let mut first: Option<(Value, Value)> = None;
    let mut identical = true;
    let mut all_verified = true;
    let mut agree = true;
    let mut last: Option<SingleRep> = None;
    let mut planned = args.repeats.max(1);
    let mut i = 0usize;
    while i < planned {
        let before = bench.measure(inst);
        let rep = single_rep(p, c, target, shape, strict)?;
        let after = bench.measure(inst);
        let unit = (before + after) / 2.0;
        let online = rep
            .pass
            .online
            .as_ref()
            .ok_or("the pass has no online interval")?;
        let setup_ns: u128 = SINGLE_SETUP_PHASES
            .iter()
            .map(|ph| rep.pass.ns.get(ph).copied().unwrap_or(0))
            .sum();
        let pass_ns: u128 = rep.pass.ns.values().sum();
        if i == 0 && pass_ns + rep.rho_setup_ns + u128::from(rep.rho.wall_ns) < FAST_NS {
            planned = args.repeats_fast.max(planned);
        }
        let rho_counts = json!({
            "steps": rep.rho.result.steps,
            "walks": rep.rho.result.walks,
            "distinguished_points": rep.rho.result.distinguished_points,
            "target_ops": rep.rho.result.target_ops,
            "recovered": rep.rho.result.recovered,
        });
        let ic_counts = json!({"pass": rep.pass.counts, "recovered": rep.pass.recovered});
        match &first {
            None => first = Some((ic_counts, rho_counts)),
            Some((ic0, rho0)) => identical &= *ic0 == ic_counts && *rho0 == rho_counts,
        }
        all_verified &= rep.pass.all_verified && online.verified && rep.rho.verified;
        agree &= online.recovered.as_ref().and_then(|d| d.to_u64()) == rep.rho.result.recovered;
        let units = |ns: f64| ns / unit;
        say(format!(
            "rep {i}: set-up {:.3} s, index calculus online {:.3} ms ({:.0} units), rho online {:.3} ms ({:.0} units), unit {unit:.2} ns",
            setup_ns as f64 / 1e9,
            online.wall_ns as f64 / 1e6,
            units(online.wall_ns as f64),
            rep.rho.wall_ns as f64 / 1e6,
            units(rep.rho.wall_ns as f64),
        ));
        reps.push(json!({
            "unit_ns_before": before, "unit_ns_after": after, "unit_ns": unit,
            "setup_phases_ns": SINGLE_SETUP_PHASES.iter().map(|ph| (ph.to_string(), rep.pass.ns.get(ph).copied().unwrap_or(0))).collect::<BTreeMap<_, _>>(),
            "setup_ns": setup_ns,
            "setup_units": units(setup_ns as f64),
            "ic_online": {
                "wall_ns": online.wall_ns,
                "units": units(online.wall_ns as f64),
                "phases_ns": online.phases_ns,
                "not_entered": online.phases_ns.iter().filter(|(_, v)| v.is_none()).map(|(k, _)| *k).collect::<Vec<_>>(),
                "outer_ns": online.outer_ns,
                "clock_lap_ns": rep.pass.ns.get("online").copied().unwrap_or(0),
                "trials": online.trials,
                "recovered": online.recovered.as_ref().map(ToString::to_string),
                "verified": online.verified,
            },
            "ic_replay_ns": online.replay_ns,
            "rho_setup_ns": rep.rho_setup_ns,
            "rho_setup_units": units(rep.rho_setup_ns as f64),
            "rho_online": {
                "wall_ns": rep.rho.wall_ns,
                "units": units(rep.rho.wall_ns as f64),
                "phases_ns": rep.rho.phases_ns,
                "outer_ns": rep.rho.outer_ns,
                "steps": rep.rho.result.steps,
                "walks": rep.rho.result.walks,
                "distinguished_points": rep.rho.result.distinguished_points,
                "rounds": rep.rho.result.rounds,
                "restart_batches": rep.rho.result.restart_batches,
                "target_ops": rep.rho.result.target_ops,
                "units_per_step": units(rep.rho.wall_ns as f64) / rep.rho.result.steps.max(1) as f64,
                "recovered": rep.rho.result.recovered,
                "verified": rep.rho.verified,
            },
            "rho_replay_ns": rep.rho.replay_ns,
            "online_speedup": rep.rho.wall_ns as f64 / online.wall_ns.max(1) as f64,
        }));
        last = Some(rep);
        i += 1;
    }
    let last = last.expect("at least one repetition");
    // The thread's model of a rho step (ledger §19–§22): one batched
    // addition and the table canonicalisation, priced in this process.
    let prices = step_prices(inst, &mut bench, args.price_rounds)?;
    let canonical_step = prices["canonical_step_units"].as_f64().unwrap_or(f64::NAN);
    let walk_ops = last
        .rho
        .result
        .counters
        .get("walk_operations")
        .copied()
        .unwrap_or(0) as f64;
    let col = |f: &dyn Fn(&Value) -> f64| -> Vec<f64> { reps.iter().map(f).collect() };
    let ic_wall = col(&|r| r["ic_online"]["wall_ns"].as_f64().unwrap_or(f64::NAN));
    let rho_wall = col(&|r| r["rho_online"]["wall_ns"].as_f64().unwrap_or(f64::NAN));
    let (ic_at, rho_at) = (median_index(&ic_wall), median_index(&rho_wall));
    let ic_rep = &reps[ic_at];
    let rho_rep = &reps[rho_at];
    let med = |f: &dyn Fn(&Value) -> f64| median(&col(f));
    let ic_units = med(&|r| r["ic_online"]["units"].as_f64().unwrap_or(f64::NAN));
    let rho_units = med(&|r| r["rho_online"]["units"].as_f64().unwrap_or(f64::NAN));
    let setup_units = med(&|r| r["setup_units"].as_f64().unwrap_or(f64::NAN));
    let rho_setup_units = med(&|r| r["rho_setup_units"].as_f64().unwrap_or(f64::NAN));
    let spread = |v: &[f64]| {
        v.iter().cloned().fold(f64::MIN, f64::max) / v.iter().cloned().fold(f64::MAX, f64::min)
    };

    // Certificates, from the last repetition: every repetition's answer
    // is the same (`counts_identical_across_repetitions`).
    let ic_d = last.pass.online.as_ref().and_then(|o| o.recovered.clone());
    let rho_d = last.rho.result.recovered.map(BigUint::from);
    let ic_method =
        "Koblitz index calculus: `ic price` single-target descent on the parameter file's recipe";
    let rho_method = last.rho.result.method.clone();
    let ic_certificate = ic_d
        .as_ref()
        .map(|d| certificate(c, "ic", ic_method, &target.q, d));
    let rho_certificate = rho_d
        .as_ref()
        .map(|d| certificate(c, "rho", &rho_method, &target.q, d));
    let target_record = json!({"curve": curve_binding(c), "point": point_pair(&target.q)});
    let columns = last.pass.counts["select"]["columns"].clone();
    let status = if all_verified && identical && agree {
        "complete"
    } else {
        "failed"
    };
    let ms = |ns: &Value| ns.as_f64().map(|v| v / 1e6);
    let target_json = json!({
        "record": target.record,
        "point": point_pair(&target.q),
        "sha256": canonical_sha256(&target_record),
        "hashed": target_record,
    });
    let recipe = json!({
        "points_requested": match &p.factor_base { FactorBaseSource::Spec { spec } => serde_json::to_value(spec).unwrap_or(Value::Null), _ => Value::Null },
        "collection_window": p.collection_window, "collection_aim": p.collection_aim,
        "unit_trials": p.collection.unit_trials, "units": p.collection.units, "max_units": p.collection.max_units,
        "summands": p.summands, "descent_summands": p.descent_summands, "max_trials": p.max_trials, "seed": p.seed,
    });
    // The median repetition's phases under the claim schema's names; a
    // phase the descent never entered is 0 and listed as such.
    let ic_phase_ms: BTreeMap<&str, f64> = IC_ONLINE_PHASES
        .iter()
        .zip(IC_CLAIM_PHASES)
        .map(|(session, (_, field))| {
            (
                field,
                ic_rep["ic_online"]["phases_ns"][*session]
                    .as_f64()
                    .unwrap_or(0.0)
                    / 1e6,
            )
        })
        .collect();
    let median_json = json!({
            "ic_online_repetition": ic_at,
            "rho_online_repetition": rho_at,
            "ic_online_wall_ns": ic_rep["ic_online"]["wall_ns"],
            "ic_online_wall_ms": ms(&ic_rep["ic_online"]["wall_ns"]),
            "ic_online_phases_ns": ic_rep["ic_online"]["phases_ns"],
            "ic_online_phases_not_entered": ic_rep["ic_online"]["not_entered"],
            "ic_online_phase_ms": ic_phase_ms,
            "rho_online_wall_ns": rho_rep["rho_online"]["wall_ns"],
            "rho_online_wall_ms": ms(&rho_rep["rho_online"]["wall_ns"]),
            "rho_online_phases_ns": rho_rep["rho_online"]["phases_ns"],
            "online_speedup": rho_rep["rho_online"]["wall_ns"].as_f64().unwrap_or(f64::NAN) / ic_rep["ic_online"]["wall_ns"].as_f64().unwrap_or(f64::NAN),
            "ic_online_units": ic_units,
            "rho_online_units": rho_units,
            "setup_units": setup_units,
            "rho_setup_units": rho_setup_units,
            "s_ic_online": ic_units / sqrt_r,
            "s_rho_online": rho_units / sqrt_r,
            "s_setup": setup_units / sqrt_r,
            "s_rho_setup": rho_setup_units / sqrt_r,
            "s_ic_cold": (setup_units + ic_units) / sqrt_r,
            "s_rho_cold": (rho_setup_units + rho_units) / sqrt_r,
            "cold_ratio_ic_over_rho": (setup_units + ic_units) / (rho_setup_units + rho_units),
            "rho_units_per_step": med(&|r| r["rho_online"]["units_per_step"].as_f64().unwrap_or(f64::NAN)),
            "rho_model_units": walk_ops * canonical_step,
            "s_rho_model": walk_ops * canonical_step / sqrt_r,
            "online_speedup_rho_model": walk_ops * canonical_step / ic_units,
            "rho_step_over_model": med(&|r| r["rho_online"]["units_per_step"].as_f64().unwrap_or(f64::NAN)) / canonical_step,
            "ic_replay_ns": med(&|r| r["ic_replay_ns"].as_f64().unwrap_or(f64::NAN)),
            "rho_replay_ns": med(&|r| r["rho_replay_ns"].as_f64().unwrap_or(f64::NAN)),
    });
    let online_interval = json!({
        "ic_start_event": "first target-dependent operation after reusable setup: the descent's first query",
        "ic_stop_event": "the descent returns a logarithm it has checked as [d]G = Q on the single-word ladder",
        "rho_start_event": "first target-dependent operation after reusable setup: lifting the target for the walk",
        "rho_stop_event": "the walk returns a logarithm it has checked as [d]G = Q on the single-word ladder",
        "ic_included_stages": IC_CLAIM_PHASES.iter().map(|(stage, _)| *stage).collect::<Vec<_>>(),
        "ic_phase_map": IC_CLAIM_PHASES.iter().zip(IC_ONLINE_PHASES).map(|((stage, _), session)| (stage.to_string(), session)).collect::<BTreeMap<_, _>>(),
        "rho_included_stages": ["walk", "collision", "recovery_check"],
        "rho_phase_map": {"rho_solve": ["walk", "collision"], "recovery_check": ["recovery_check"]},
        "replay": "[d]G = Q in the general arithmetic, for both arms after their intervals; timed as ic_replay_ns and rho_replay_ns",
    });
    let rho_policy = json!({
        "worker_count": threads,
        "lanes": lanes,
        "distinguished_point_bits": dp_bits,
        "seed": seed,
        "walk_policy": rho_method,
        "collision_policy": "distinguished points in one table for the target's walks; a stored point reached with a different coefficient of Q gives d, checked as [d]G = Q; fruitless short cycles are left by doubling their least point",
        "distinguished_point_memory_bytes": last.rho.result.counters.get("distinguished_point_table_bytes"),
        "counters": last.rho.result.counters,
    });
    let resource_envelope = json!({
        "rayon_threads": threads,
        "cpus_allowed_list": cpus_allowed_list(),
        "arms": "one process, both arms rebuilt from nothing each repetition, the index calculus first",
        "os": std::env::consts::OS,
        "arch": std::env::consts::ARCH,
    });
    let certificates = json!({
        "ic": ic_certificate,
        "ic_sha256": ic_certificate.as_ref().map(canonical_sha256),
        "rho": rho_certificate,
        "rho_sha256": rho_certificate.as_ref().map(canonical_sha256),
    });
    let identity_inputs = json!({
        "fixture": fixture(c, target),
        "factor_base_orbits": last.pass.factor_base_orbits,
        "columns": columns,
        "column_convention": "representative",
    });
    Ok(json!({
        "schema_version": 1,
        "operation": "price_single_target",
        "status": status,
        "params": args.params.display().to_string(),
        "name": p.name,
        "curve": inst.name, "n": inst.n, "a": p.curve.curve_a, "r": r, "log2_r": (r as f64).log2(),
        "target_count": 1,
        "target": target_json,
        "recipe": recipe,
        "threads": threads,
        "single_thread": threads == 1,
        "strict_phase_session": strict,
        "counts": first.as_ref().map(|f| f.0.clone()),
        "rho_counts": first.as_ref().map(|f| f.1.clone()),
        "counts_identical_across_repetitions": identical,
        "all_verified": all_verified,
        "ic_and_rho_agree": agree,
        "repetitions": reps,
        "median": median_json,
        "spread_max_over_min": {
            "ic_online": spread(&ic_wall),
            "rho_online": spread(&rho_wall),
        },
        "online_interval": online_interval,
        "step_prices": prices,
        "rho_model": "rho's counted walk operations at the canonical step, one batched addition plus the table canonicalisation, priced in this process (ledger §19–§22); a model beside the measured interval, not a measurement of it",
        "rho_setup": {
            "what": "the curve's single-word arithmetic, the signed-Frobenius classes and the jump table [c_j]G: nothing depends on the target",
            "setup_ops": last.rho.result.setup_ops,
        },
        "rho_policy": rho_policy,
        "resource_envelope": resource_envelope,
        "certificates": certificates,
        "identity_inputs": identity_inputs,
        "commit": super::git_commit(),
        "binary_blake3": super::binary_hash(),
    }))
}

/// Canonical-step and Bailey-step prices, in the unit, measured in
/// interleaved rounds (ledger §19.1 b).
fn step_prices(
    inst: &BinaryInstance,
    bench: &mut UnitBench,
    rounds: usize,
) -> Result<Value, String> {
    let classes = SignedFrobeniusClasses::new(inst).ok_or("no signed-Frobenius classes")?;
    let canon = FrobeniusCanon::new(&inst.fast.field, inst.n).ok_or("no normal basis")?;
    let powers = FrobeniusPowers::new(&inst.fast.field, inst.n);
    let bg = BinaryGroup(&inst.fast);
    let pts = bench.batch.clone();
    let mut sink = 0u64;
    let (mut canon_units, mut bailey_units, mut units) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..rounds {
        let unit = bench.measure(inst);
        units.push(unit);
        let t = Instant::now();
        for i in 0..CANONS {
            let (rep, mu) = classes.canon(&bg, pts[i & 1023]);
            sink ^= rep.x ^ rep.y ^ mu;
        }
        canon_units.push(t.elapsed().as_nanos() as f64 / CANONS as f64 / unit);
        let t = Instant::now();
        for i in 0..CANONS {
            let p = pts[i & 1023];
            let j = 3 + (canon.coords(p.x).count_ones() / 2) % 8;
            sink ^= powers.apply(j, p.x) ^ powers.apply(j, p.y);
        }
        bailey_units.push(t.elapsed().as_nanos() as f64 / CANONS as f64 / unit);
    }
    std::hint::black_box(sink);
    let canon_m = median(&canon_units);
    let bailey_m = median(&bailey_units);
    Ok(json!({
        "rounds": rounds,
        "unit_ns": units,
        "canonicalisation_units": canon_units,
        "bailey_overhead_units": bailey_units,
        "canonical_step_units": 1.0 + canon_m,
        "bailey_step_units": 1.0 + bailey_m,
    }))
}

fn rho_json(batch: &RhoBatchResult, r: u64) -> Value {
    json!({
        "targets": batch.targets,
        "gae": batch.gae,
        "setup_gae": batch.setup_ops.gae(),
        "s_counted_per_target": batch.gae / (batch.targets.max(1) as f64 * (r as f64).sqrt()),
        "all_verified": batch.all_verified,
        "solved_by": batch.per_target.iter().map(|t| t.solved_by).collect::<Vec<_>>(),
        "steps_per_target": batch.per_target.iter().map(|t| t.steps).collect::<Vec<_>>(),
        "recovered": batch.per_target.iter().map(|t| t.recovered).collect::<Vec<_>>(),
        "counters": batch.counters,
        "wall_ns": batch.wall_ns,
    })
}

pub fn run(args: PriceArgs, quiet: bool) -> Result<Value, String> {
    let threads = rayon::current_num_threads();
    if threads != 1 && !args.allow_threads {
        return Err(format!(
            "ic price times phases on one thread; set RAYON_NUM_THREADS=1 (found {threads}) or pass --allow-threads"
        ));
    }
    let say = |msg: String| {
        if !quiet {
            eprintln!("{msg}");
        }
    };
    let p = workflow::load_params(&args.params)?;
    if p.curve.subfield != 1 || p.curve.curve_b != 1 || p.curve.curve_a > 1 {
        return Err("ic price prices Koblitz curves (subfield 1, b = 1, a in {0, 1})".into());
    }
    let inst = koblitz_instance(p.curve.curve_a, p.curve.degree)
        .ok_or("no single-word Koblitz instance for this curve")?;
    let c = experiment::curve(
        p.curve.degree,
        p.curve.curve_a,
        p.curve.subfield,
        p.curve.curve_b,
    )?;
    if c.subgroup_order.to_u64() != Some(inst.r) {
        return Err("the workflow's curve and the rho instance disagree on r".into());
    }
    let r = inst.r;
    let sqrt_r = (r as f64).sqrt();
    let k = p.targets.len();
    if k == 0 {
        return Err("no targets".into());
    }
    // Targets, outside every clock.
    let targets: Vec<Target> = p
        .targets
        .iter()
        .map(|t| {
            workflow::resolve_target(&c, t).map(|(q, expected, record)| Target {
                q,
                expected,
                record,
            })
        })
        .collect::<Result<_, _>>()?;
    if args.single_target {
        if k != 1 {
            return Err(format!(
                "--single-target prices exactly one target; the file names {k}"
            ));
        }
        return run_single(&args, &p, &c, &inst, &targets[0], threads, true, &say);
    }
    let mut bench = UnitBench::new(&inst);
    // One untimed pass of the unit so its own first touch is not a price.
    bench.measure(&inst);

    let mut reps: Vec<Value> = Vec::new();
    let mut first_counts: Option<Value> = None;
    let mut recovered_first: Option<Vec<Option<String>>> = None;
    let mut counts_identical = true;
    let mut all_verified = true;
    let mut planned = args.repeats.max(1);
    let mut i = 0usize;
    while i < planned {
        let before = bench.measure(&inst);
        let pass = one_pass(&p, &targets)?;
        let after = bench.measure(&inst);
        let unit = (before + after) / 2.0;
        let total_ns: u128 = pass.ns.values().sum();
        if i == 0 && total_ns < FAST_NS {
            planned = args.repeats_fast.max(planned);
        }
        match &first_counts {
            None => {
                first_counts = Some(pass.counts.clone());
                recovered_first = Some(pass.recovered.clone());
            }
            Some(first) => {
                counts_identical &=
                    *first == pass.counts && recovered_first.as_ref() == Some(&pass.recovered);
            }
        }
        all_verified &= pass.all_verified;
        let units: BTreeMap<&str, f64> = PHASES
            .iter()
            .map(|ph| (*ph, pass.ns.get(ph).copied().unwrap_or(0) as f64 / unit))
            .collect();
        let total_units = total_ns as f64 / unit;
        say(format!(
            "rep {i}: {:.3} s, unit {unit:.2} ns, {total_units:.0} units, S {:.4} per target",
            total_ns as f64 / 1e9,
            total_units / (k as f64 * sqrt_r)
        ));
        reps.push(json!({
            "unit_ns_before": before, "unit_ns_after": after, "unit_ns": unit,
            "phases_ns": PHASES.iter().map(|ph| (ph.to_string(), pass.ns.get(ph).copied().unwrap_or(0))).collect::<BTreeMap<_, _>>(),
            "phases_units": units,
            "total_ns": total_ns,
            "total_units": total_units,
            "s_per_target": total_units / (k as f64 * sqrt_r),
        }));
        i += 1;
    }
    let col = |key: &str| -> Vec<f64> {
        reps.iter()
            .map(|r| r[key].as_f64().unwrap_or(f64::NAN))
            .collect()
    };
    let totals = col("total_units");
    let spread = totals.iter().cloned().fold(f64::MIN, f64::max)
        / totals.iter().cloned().fold(f64::MAX, f64::min);
    let median_phases: BTreeMap<&str, f64> = PHASES
        .iter()
        .map(|ph| {
            let v: Vec<f64> = reps
                .iter()
                .map(|r| r["phases_units"][*ph].as_f64().unwrap_or(0.0))
                .collect();
            (*ph, median(&v))
        })
        .collect();
    let median_total = median(&totals);
    let s_ic = median_total / (k as f64 * sqrt_r);
    let counts = first_counts.unwrap_or(Value::Null);
    // Cold: the shared phases and one mean descent, from this run.
    let per_target_phases = ["descent", "verify_final"];
    let shared: f64 = median_phases
        .iter()
        .filter(|(ph, _)| !per_target_phases.contains(ph))
        .map(|(_, v)| v)
        .sum();
    let per_target: f64 = per_target_phases
        .iter()
        .map(|ph| median_phases[ph])
        .sum::<f64>()
        / k as f64;
    let s_cold = (shared + per_target) / sqrt_r;
    // Read-outs of the same phases, never a different total: the work a
    // model of the method prices, the constructions around it, and the
    // checks.
    let group = |names: &[&str]| -> f64 { names.iter().map(|ph| median_phases[ph]).sum() };
    let groups = json!({
        "work": group(&["select", "build", "collect", "la", "descent"]),
        "constructions": group(&["setup", "select_projection", "collect_setup", "logs_setup", "descent_setup", "other"]),
        "verification": group(&["verify", "verify_final"]),
        "what": "work = the phases the frozen model prices; constructions = curve, projections, collectors, coverage and solvers built around them; verification = relation and final checks. The three sum to the median phases, which need not equal the median total.",
    });

    let floor_k = super::rho::batch_law(k) * (std::f64::consts::PI / (4.0 * inst.n as f64)).sqrt();
    let mut report = json!({
        "schema_version": 1,
        "operation": "price",
        "status": if all_verified && counts_identical { "complete" } else { "failed" },
        "params": args.params.display().to_string(),
        "name": p.name,
        "curve": inst.name, "n": inst.n, "a": p.curve.curve_a, "r": r, "log2_r": (r as f64).log2(),
        "targets": k,
        "recipe": {
            "points_requested": match &p.factor_base { FactorBaseSource::Spec { spec } => serde_json::to_value(spec).unwrap_or(Value::Null), _ => Value::Null },
            "collection_window": p.collection_window, "collection_aim": p.collection_aim,
            "unit_trials": p.collection.unit_trials, "units": p.collection.units, "max_units": p.collection.max_units,
            "summands": p.summands, "descent_summands": p.descent_summands, "max_trials": p.max_trials, "seed": p.seed,
        },
        "threads": threads,
        "single_thread": threads == 1,
        "counts": counts,
        "counts_identical_across_repetitions": counts_identical,
        "all_verified": all_verified,
        "recovered": recovered_first,
        "repetitions": reps,
        "first_repetition_total_units": totals.first(),
        "spread_max_over_min": spread,
        "median": {
            "phases_units": median_phases,
            "total_units": median_total,
            "s_per_target": s_ic,
            "s_cold": s_cold,
            "shared_units": shared,
            "per_target_units": per_target,
            "groups": groups,
        },
        "floor": {"per_target_k": floor_k, "batch_law": super::rho::batch_law(k)},
        "commit": super::git_commit(),
    });

    if let Some(seed) = args.rho_seed {
        let prices = step_prices(&inst, &mut bench, args.price_rounds)?;
        let step = prices["canonical_step_units"].as_f64().unwrap_or(f64::NAN);
        let bailey = prices["bailey_step_units"].as_f64().unwrap_or(f64::NAN);
        let q_fast: Vec<FastPoint> = targets.iter().map(|t| inst.fast.lift(&t.q)).collect();
        say(format!("batch rho on the same {k} targets …"));
        let batch = signed_frobenius_rho_batch(&inst, &q_fast, seed)
            .ok_or("no signed-Frobenius classes")?;
        // The same answers as the index calculus, target by target.
        let agree =
            batch
                .per_target
                .iter()
                .zip(&targets)
                .all(|(b, t)| match (&t.expected, b.recovered) {
                    (Some(e), Some(d)) => e.to_u64() == Some(d),
                    (None, Some(_)) => b.verified,
                    _ => false,
                });
        let s_counted = batch.gae / (k as f64 * sqrt_r);
        let mut rho = rho_json(&batch, r);
        rho["seed"] = json!(seed);
        rho["agrees_with_expected"] = json!(agree);
        rho["s_priced_per_target"] = json!(s_counted * step);
        rho["s_bailey_model_per_target"] = json!(s_counted * bailey);
        report["rho_batch"] = rho;
        report["step_prices"] = prices;
        report["ratios"] = json!({
            "ic_over_batch_rho": s_ic / (s_counted * step),
            "ic_over_batch_rho_bailey_model": s_ic / (s_counted * bailey),
            "ic_over_floor": s_ic / floor_k,
            "batch_rho_counted_over_floor": s_counted / floor_k,
        });
        if !(batch.all_verified && agree) {
            report["status"] = json!("failed");
        }
        if let Some(cold_seed) = args.cold_rho_seed {
            say(format!("single-target rho on each of the {k} targets …"));
            let mut gae = Vec::with_capacity(k);
            let mut ok = true;
            for (i, t) in q_fast.iter().enumerate() {
                let one = signed_frobenius_rho_batch(
                    &inst,
                    std::slice::from_ref(t),
                    cold_seed + i as u64,
                )
                .ok_or("no signed-Frobenius classes")?;
                ok &= one.all_verified
                    && targets[i].expected.as_ref().and_then(|e| e.to_u64())
                        == one.per_target[0].recovered;
                gae.push(one.gae);
            }
            let mean = gae.iter().sum::<f64>() / k as f64;
            let sd = (gae.iter().map(|g| (g - mean).powi(2)).sum::<f64>() / (k.max(2) - 1) as f64)
                .sqrt();
            let s_cold_counted = mean / sqrt_r;
            report["rho_cold"] = json!({
                "seed_base": cold_seed,
                "gae_per_target": gae,
                "s_counted": s_cold_counted,
                "s_counted_sd_of_mean": sd / (k as f64).sqrt() / sqrt_r,
                "s_priced": s_cold_counted * step,
                "all_verified": ok,
            });
            report["ratios"]["cold_ic_over_single_rho"] = json!(s_cold / (s_cold_counted * step));
            if !ok {
                report["status"] = json!("failed");
            }
        }
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(targets: usize, descent: u8) -> workflow::WorkflowParams {
        serde_json::from_value(json!({
            "schema_version": 1,
            "curve": {"degree": 17, "curve_a": 1},
            "summands": 3,
            "descent_summands": descent,
            "collection_window": 8,
            "collection_aim": true,
            "seed": 5,
            "max_trials": 100000,
            "collection": {"unit_trials": 16, "units": 1, "max_units": 1000},
            "factor_base": {"mode": "spec", "spec": {"kind": "subgroup_orbits", "points": 272, "seed": 5}},
            "targets": (0..targets).map(|i| json!({"random_seed": 500 + i})).collect::<Vec<_>>(),
        }))
        .expect("parameters parse")
    }

    fn targets(p: &workflow::WorkflowParams) -> Vec<Target> {
        let c = experiment::curve(17, 1, 1, 1).expect("curve");
        p.targets
            .iter()
            .map(|t| {
                let (q, expected, record) = workflow::resolve_target(&c, t).expect("target");
                Target {
                    q,
                    expected,
                    record,
                }
            })
            .collect()
    }

    fn one_public_target(descent: u8) -> (workflow::WorkflowParams, Target) {
        let mut p = params(1, descent);
        p.targets = serde_json::from_value(json!([{"public_hash_seed": 2026093001u64}]))
            .expect("a public target");
        let c = experiment::curve(17, 1, 1, 1).expect("curve");
        let (q, expected, record) = workflow::resolve_target(&c, &p.targets[0]).expect("target");
        assert!(expected.is_none() && !record.target_scalar_constructed);
        (
            p,
            Target {
                q,
                expected,
                record,
            },
        )
    }

    #[test]
    fn a_single_target_pass_splits_its_online_interval_exactly() {
        let c = experiment::curve(17, 1, 1, 1).expect("curve");
        let inst = koblitz_instance(1, 17).expect("a Koblitz instance");
        for descent in [2u8, 3] {
            let (p, target) = one_public_target(descent);
            let lanes = ParallelRho::default_lanes(&inst);
            let shape = RhoShape {
                seed: 99,
                lanes,
                dp_bits: ParallelRho::default_dp_bits(&inst, lanes),
            };
            let first = single_rep(&p, &c, &target, shape, false).expect("a repetition");
            let second = single_rep(&p, &c, &target, shape, false).expect("a repetition");
            for rep in [&first, &second] {
                let online = rep.pass.online.as_ref().expect("an online interval");
                assert!(online.verified && rep.rho.verified, "m = {descent}");
                // Every phase of the descent was entered, and together
                // they are the interval, to the nanosecond.
                assert!(
                    online.phases_ns.values().all(Option::is_some),
                    "{:?}",
                    online.phases_ns
                );
                assert_eq!(
                    online.phases_ns.values().flatten().sum::<u64>(),
                    online.wall_ns
                );
                assert!(online.wall_ns > 0 && u128::from(online.wall_ns) <= online.outer_ns);
                assert_eq!(
                    rep.rho.phases_ns.values().flatten().sum::<u64>(),
                    rep.rho.wall_ns
                );
                assert!(rep.rho.phases_ns["rho_solve"].is_some());
                // Both arms found the same logarithm of the same point.
                assert_eq!(
                    online.recovered.as_ref().and_then(|d| d.to_u64()),
                    rep.rho.result.recovered
                );
                // The set-up is on the exclusive clock, apart from the target.
                for phase in [
                    "select",
                    "build",
                    "collect",
                    "la",
                    "descent_setup",
                    "online",
                    "replay",
                ] {
                    assert!(rep.pass.ns.contains_key(phase), "{phase} was not clocked");
                }
                assert!(!rep.pass.ns.contains_key("descent"));
                let orbits = rep
                    .pass
                    .factor_base_orbits
                    .as_ref()
                    .expect("the base's orbits");
                assert_eq!(
                    json!(orbits.len()),
                    rep.pass.counts["select"]["signed_orbits"]
                );
            }
            assert_eq!(first.pass.counts, second.pass.counts);
            assert_eq!(first.rho.result.steps, second.rho.result.steps);
            assert_eq!(first.rho.result.target_ops, second.rho.result.target_ops);
        }
    }

    #[test]
    fn certificates_hash_the_canonical_form() {
        // serde_json's compact form with sorted keys is the canonical
        // form the replay recomputes; a reordered map must hash the same.
        let a = json!({"b": 1, "a": [2, "x"]});
        assert_eq!(serde_json::to_string(&a).unwrap(), r#"{"a":[2,"x"],"b":1}"#);
        assert_eq!(
            canonical_sha256(&a),
            hex::encode(sha256(br#"{"a":[2,"x"],"b":1}"#))
        );
        let c = experiment::curve(17, 1, 1, 1).expect("curve");
        let (_, target) = one_public_target(2);
        let d = BigUint::from(12345u32);
        let one = certificate(&c, "ic", "m", &target.q, &d);
        assert_eq!(canonical_sha256(&one), canonical_sha256(&one.clone()));
        assert_ne!(
            canonical_sha256(&one),
            canonical_sha256(&certificate(&c, "rho", "m", &target.q, &d))
        );
    }

    #[test]
    fn a_pass_solves_every_target_and_repeats_its_counts() {
        for descent in [2u8, 3] {
            let p = params(4, descent);
            let t = targets(&p);
            let first = one_pass(&p, &t).expect("a pass");
            let second = one_pass(&p, &t).expect("a pass");
            assert!(first.all_verified && second.all_verified, "m = {descent}");
            assert_eq!(first.counts, second.counts, "m = {descent}");
            assert_eq!(first.recovered, second.recovered, "m = {descent}");
            assert!(first.recovered.iter().all(Option::is_some));
            // Every clock is a declared phase, and the work is on its own.
            assert!(first.ns.keys().all(|k| PHASES.contains(k)));
            for phase in ["select", "build", "collect", "la", "descent"] {
                assert!(first.ns.contains_key(phase), "{phase} was not clocked");
            }
            assert_eq!(first.counts["descent"]["summands"], json!(descent));
        }
    }
}
