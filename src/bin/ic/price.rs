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
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use clap::Args;
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_instance, signed_frobenius_rho_batch, BinaryGroup, BinaryInstance, CountedGroup,
    GroupOps, RhoBatchResult, RhoClasses, SignedFrobeniusClasses,
};
use crypto_lib::cryptanalysis::koblitz_fast::{
    BatchScratch, FastPoint, FrobeniusCanon, FrobeniusPowers,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    projected_signed_orbit_count, CollectedRelation, ColumnCoverage, DecompositionStrategy,
    FactorBaseLogSolver, IndividualLogSolver, RelationCollector, RelationWorkUnit,
};
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

/// One target, as the method receives it, with the answer it must give.
struct Target {
    q: BinaryPoint,
    expected: Option<BigUint>,
}

/// One pass of the pipeline: its phase clocks and its native counts.
struct Pass {
    ns: BTreeMap<&'static str, u128>,
    counts: Value,
    recovered: Vec<Option<String>>,
    all_verified: bool,
}

fn one_pass(p: &workflow::WorkflowParams, targets: &[Target]) -> Result<Pass, String> {
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
    let mut trials = Vec::with_capacity(targets.len());
    let mut recovered = Vec::with_capacity(targets.len());
    let mut all_verified = true;
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
    })
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
        .map(|t| workflow::resolve_target(&c, t).map(|(q, expected, _)| Target { q, expected }))
        .collect::<Result<_, _>>()?;
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
                let (q, expected, _) = workflow::resolve_target(&c, t).expect("target");
                Target { q, expected }
            })
            .collect()
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
