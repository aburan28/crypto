//! `pdp-split`: the split and hybrid point-decomposition systems of
//! `crypto_lib::cryptanalysis::pdp_split` as an `ic.pipeline` oracle
//! (`research/pdp_split_slots_20261012/PROTOCOL.md`).
//!
//! ```text
//! pdp-split:m=M,plan=P,guess=G,slots=one|all,max_free=F
//! ```
//!
//! - `plan` is the block plan joined by `+` (`3`, `2+1`, `2+2`, `2+1+1`);
//!   its sum must be `m`.
//! - `guess` fixes the first `G` Boolean unknowns to each of `2^G` values
//!   (Bettale–Faugère–Perret), branch by branch, until one lifts.
//! - `slots=one` solves the one system with every summand over `V`;
//!   `slots=all` (on `koblitz-closure`) solves, in a fixed order, one system
//!   per multiset of rotations `(0, j₂ ≤ … ≤ j_m)` with summand `i` over
//!   `σ^{jᵢ}(V)`, until one lifts.
//! - `max_free` caps the free unknowns the read-out enumerates.
//!
//! Every parameter that changes the search is required.  F4 runs with no
//! wall-clock budget, so a query is a deterministic function of its target.
//!
//! Charged: the group additions of the sign lift
//! (`ic_boundary::lift_abscissae`).  Counted, not charged: F4's word XORs,
//! under the `f4-f2` adapter's unit (elimination only, so not a whole-run
//! count and never priced), with branches, systems, solutions and the
//! largest solving degree in `solver.extra`.  `S` is a lower bound.  The
//! whole query — the symbolic build of each split system, each branch's
//! specialisation, F4 end to end and the read-out — is also reported in
//! `bitop/v1` primitives through the `bitop_*` keys of `solver.extra`
//! with `bitop_complete = 1` (`crypto_lib::cryptanalysis::bitop_tally`).
//!
//! The `koblitz-closure` base (`ic_framework::plugins::KoblitzClosureBase`)
//! is resolved here too, so `ic.pipeline` can reach it from ecbench.

use std::collections::BTreeMap;
use std::time::Instant;

use crypto_lib::cryptanalysis::bitop_tally;
use crypto_lib::cryptanalysis::ic_boundary::{
    lift_abscissae, BinaryGroup, BinaryInstance, FactorBase, GroupOps, OracleCounters,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::KoblitzClosureBase;
use crypto_lib::cryptanalysis::ic_framework::stages::{
    DecompositionOracle, FactorBaseBuilder, InstanceCtx, Params, SolverCost, SolverTotals,
    SystemShape,
};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::pdp_split::{
    build_split_system_bases, lift_x, plan_label, rotate_words, slot_tuples, solve_branch,
    span_words, summands_of_bases,
};

/// The `f4-f2` adapter's unit: F4's elimination XORs, not a whole-run count.
const F4_UNIT: &str = "word XORs (elimination only)";

/// Environment prefixes the F4 code reads for tuning overrides.
const TUNING_PREFIXES: [&str; 3] = ["KIC_", "F4_", "SOLVER_"];

pub struct PdpSplitOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    plan: Vec<usize>,
    guess: usize,
    max_free: usize,
    /// Rotation tuples tried per query, in order.
    slots: Vec<Vec<u32>>,
    /// `V`'s basis words.
    v: Vec<u64>,
    st: Option<FieldStructure>,
    shape: Option<SystemShape>,
    totals: SolverTotals,
}

impl<'i> PdpSplitOracle<'i> {
    pub fn new(summands: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            summands,
            instance,
            plan: Vec::new(),
            guess: 0,
            max_free: 0,
            slots: Vec::new(),
            v: Vec::new(),
            st: None,
            shape: None,
            totals: SolverTotals::default(),
        }
    }
}

fn required<'p>(params: &'p Params, key: &str, hint: &str) -> Result<&'p str, String> {
    params
        .get(key)
        .ok_or_else(|| format!("missing parameter `{key}`; {hint}"))
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for PdpSplitOracle<'_> {
    fn name(&self) -> &str {
        "pdp-split"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "Weil-descend the {} plan of {} summands, guess {} unknowns, {} slots; F4",
            params.get("plan").unwrap_or("?"),
            self.summands,
            params.get("guess").unwrap_or("?"),
            params.get("slots").unwrap_or("?"),
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands; the plan's sum"),
            (
                "plan",
                "block plan joined by '+': 3, 2+1, 2+2, 2+1+1 (m = 2: 2)",
            ),
            (
                "guess",
                "Boolean unknowns fixed per branch, 2^guess branches",
            ),
            (
                "slots",
                "one (every summand over V) or all (every rotation multiset over koblitz-closure)",
            ),
            ("max_free", "free unknowns the solution read-out enumerates"),
            ("dimension", "the base's; read from the base"),
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
                "pdp-split refuses to run with the solver override `{k}` set: it would change \
                 the search without changing the method's identity"
            ));
        }
        let inst = self.instance;
        if inst.koblitz.is_none() || inst.b != 1 {
            return Err("pdp-split needs a Koblitz instance (b = 1, for the S4 link)".into());
        }
        let plan_s = required(params, "plan", "try 3 or 2+1")?;
        self.plan = plan_s
            .split('+')
            .map(|p| {
                p.parse::<usize>()
                    .map_err(|_| format!("plan `{plan_s}` is not sizes joined by '+'"))
            })
            .collect::<Result<_, _>>()?;
        let m: usize = self.plan.iter().sum();
        if m != self.summands as usize {
            return Err(format!(
                "plan {plan_s} has {m} summands but m = {}",
                self.summands
            ));
        }
        let g = required(params, "guess", "try 0, 2 or 4")?;
        self.guess = g
            .parse()
            .map_err(|_| format!("guess `{g}` is not a number"))?;
        let mf = required(params, "max_free", "the study ran 24")?;
        self.max_free = mf
            .parse()
            .map_err(|_| format!("max_free `{mf}` is not a number"))?;
        let basis = fb.subspace_basis.clone().ok_or(
            "pdp-split needs a subspace base (koblitz-standard-subspace or koblitz-closure)",
        )?;
        let ell = basis.len();
        if self.guess > m * ell || self.guess > 20 {
            return Err(format!(
                "guess {} exceeds the {} summand unknowns (or 20)",
                self.guess,
                m * ell
            ));
        }
        self.slots = match required(params, "slots", "one or all")? {
            "one" => vec![vec![0u32; m]],
            "all" => {
                // The rotations must be in the base, or no hit can lift.
                // Every liftable abscissa of σ(V) must be a base abscissa,
                // and some must lie outside V, or rotated hits cannot lift.
                let v_span = span_words(&basis);
                // (x = 0, the 2-torsion point, is in no closure base.)
                let liftable: Vec<u64> = span_words(&rotate_words(inst, &basis, 1))
                    .into_iter()
                    .filter(|&x| x != 0 && lift_x(inst, x).is_some())
                    .collect();
                let holds = liftable.iter().all(|x| fb.x_index.contains_key(x))
                    && liftable.iter().any(|x| !v_span.contains(x));
                if !holds {
                    return Err("slots=all needs the koblitz-closure base: the run's base \
                                does not hold sigma(V)"
                        .into());
                }
                slot_tuples(inst.n, m)
            }
            other => return Err(format!("slots must be one or all, got `{other}`")),
        };
        self.v = basis;
        let k = self.plan.len();
        let n_vars = m * ell + (k - 1) * inst.n as usize;
        if n_vars > 64 {
            return Err(format!(
                "plan {plan_s} at dimension {ell} has {n_vars} > 64 unknowns"
            ));
        }
        self.st = Some(FieldStructure::new(inst.n, &inst.irreducible));
        self.shape = Some(SystemShape {
            n_vars,
            n_equations: k * inst.n as usize,
            degrees: Vec::new(),
            semi_regular_degree: None,
        });
        self.totals = SolverTotals {
            solver: format!(
                "pdp-split {} guess {} slots {} (F4, max_free {})",
                plan_label(&self.plan),
                self.guess,
                params.get("slots").unwrap_or("?"),
                self.max_free
            ),
            ..Default::default()
        };
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
        let (Some(st), Some(shape)) = (&self.st, &self.shape) else {
            return None;
        };
        if point.infinity {
            return None;
        }
        let inst = self.instance;
        let b = inst.gf.to_element(inst.b);
        let xr = inst.gf.to_element(point.x);
        let started = Instant::now();
        let tally_before = bitop_tally::snapshot();
        let mut performed = 0u64;
        let (mut xors, mut systems, mut branches, mut solutions) = (0u64, 0u64, 0u64, 0u64);
        let (mut sd_max, mut deg_max, mut timed_out) = (0u32, 0u32, false);
        let mut found = None;
        'slots: for slots in &self.slots {
            let words: Vec<Vec<u64>> = slots
                .iter()
                .map(|&j| rotate_words(inst, &self.v, j))
                .collect();
            // `bitop/v1`: each rotation is `j` squarings a basis word, and
            // every basis word is converted to a field element
            bitop_tally::field_ops(
                0,
                0,
                slots.iter().map(|&j| u64::from(j)).sum::<u64>() * self.v.len() as u64,
                0,
            );
            bitop_tally::word_ops((slots.len() * self.v.len()) as u64);
            let bases: Vec<_> = words
                .iter()
                .map(|w| w.iter().map(|&x| inst.gf.to_element(x)).collect::<Vec<_>>())
                .collect();
            let Ok(sys) = build_split_system_bases(&self.plan, &bases, &xr, &b, st) else {
                break;
            };
            systems += 1;
            for guess in 0..(1u64 << self.guess) {
                let bo = solve_branch(&sys, self.guess, guess, None, self.max_free);
                branches += 1;
                xors += bo.word_xors;
                performed += bo.word_xors_performed;
                sd_max = sd_max.max(bo.solving_degree);
                deg_max = deg_max.max(bo.degree_reached);
                timed_out |= bo.timed_out;
                let Some(sols) = bo.solutions else {
                    counters.unliftable_systems += 1;
                    continue;
                };
                for a in sols {
                    solutions += 1;
                    // one test and at most one `F_{2^n}` addition a
                    // coordinate of every summand
                    bitop_tally::word_ops((sys.m * sys.ell) as u64);
                    bitop_tally::field_ops((sys.m * sys.ell) as u64, 0, 0, 0);
                    let xs = summands_of_bases(&sys, &words, a);
                    match lift_abscissae(ctx.group, fb, ops, &xs, point) {
                        Some(idx) => {
                            found = Some(idx);
                            break 'slots;
                        }
                        None => counters.lift_failures += 1,
                    }
                }
            }
        }
        let mut extra = BTreeMap::new();
        extra.insert("systems".to_string(), systems);
        extra.insert("branches".to_string(), branches);
        extra.insert("algebraic_solutions".to_string(), solutions);
        extra.insert("solving_degree_max".to_string(), u64::from(sd_max));
        // The whole query in `bitop/v1` primitives: the XORs F4's
        // eliminations performed (`ops` counts the logical ones), and
        // everything else `bitop_tally` saw on this thread — the symbolic
        // `F_{2^n}` build of each split system, the specialisation of each
        // branch, F4's products, symbolic preprocessing, pair selection,
        // sorts, packing, lead scans and read-back, and each branch's
        // read-out.  The sign lift's group additions are charged to `ops`
        // of the group, not here.
        bitop_tally::snapshot()
            .since(tally_before)
            .put_complete(&mut extra, performed);
        let cost = SolverCost {
            ops: xors,
            op_unit: F4_UNIT.into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: Some(deg_max),
            solving_degree: Some(sd_max),
            timed_out,
            extra,
        };
        self.totals.absorb(shape, Some(&cost), timed_out);
        found
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

/// The `koblitz-closure` base, by name.
pub fn base<'i>(
    name: &str,
    instance: &'i BinaryInstance,
) -> Option<Box<dyn FactorBaseBuilder<BinaryGroup<'i>> + 'i>> {
    (name == "koblitz-closure").then(|| {
        Box::new(KoblitzClosureBase { instance })
            as Box<dyn FactorBaseBuilder<BinaryGroup<'i>> + 'i>
    })
}

/// The `pdp-split` oracle, by name.
pub fn oracle<'i>(
    name: &str,
    summands: u32,
    instance: &'i BinaryInstance,
) -> Option<Box<dyn DecompositionOracle<BinaryGroup<'i>> + 'i>> {
    (name == "pdp-split").then(|| {
        Box::new(PdpSplitOracle::new(summands, instance))
            as Box<dyn DecompositionOracle<BinaryGroup<'i>> + 'i>
    })
}
