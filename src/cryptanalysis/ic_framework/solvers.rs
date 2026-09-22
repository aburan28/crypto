//! # Polynomial-system solvers, as plug-ins.
//!
//! The [`SystemSolver`] plug point takes a boolean system and returns
//! its solutions with a cost record.  It knows nothing about elliptic
//! curves, which is the point: a Gröbner implementer can work here
//! without reading a line of curve arithmetic, and the framework can
//! ask "what does this pipeline cost if the decomposition oracle solves
//! its systems with *that* engine instead" without touching anything
//! else.
//!
//! ## What ships here
//!
//! | name | engine | native unit |
//! |:--|:--|:--|
//! | `buchberger-f2` | boolean-ring Buchberger with the coprime and chain criteria | monomial operations |
//! | `xl-f2` | XL: multiply out to a degree, then linearise | monomial operations |
//! | `sat-cdcl` | CDCL with Tseitin monomials and native parity rows | conflicts |
//! | `exhaustive` | evaluate every equation at every point | monomial tests |
//!
//! `exhaustive` is not a strawman.  It is the **reference** the others
//! are measured against, in the sense `AGENTS.md` §1 means: the best
//! algorithm that already solves the same problem, priced in the same
//! unit on the same instances.  A solver that cannot beat it on a cell
//! has not earned its place in that cell, however good its asymptotics
//! are said to be.
//!
//! ## Adding an engine — F4, F5, or anything else
//!
//! Implement [`SystemSolver`] and register it in [`solver_registry`].
//! The contract is short and all of it matters:
//!
//! - **Count something.**  `SolverCost::ops` with an `op_unit` naming
//!   it.  Wall time is never the metric (`AGENTS.md` §6); an engine
//!   that reports only time cannot be compared across hosts and will
//!   not be accepted into the table as a speed.
//! - **Respect the budget.**  Return [`SolverVerdict::BudgetExceeded`],
//!   never a wrong answer and never an unbounded run.
//! - **Never guess `Unsatisfiable`.**  Say it only when the system is
//!   decided.  An engine that reports a timeout as "no solution"
//!   quietly lowers the hit rate of every relation phase built on it.
//! - **Report a degree if you have one.**  `solving_degree` is the
//!   highest degree at which the run learned something; it is the one
//!   comparable against a degree bound.  `degree_reached` may be higher
//!   and is a property of the strategy.
//!
//! See `docs/ic/FRAMEWORK.md` for a worked example of adding one.

use std::collections::BTreeMap;
use std::time::{Duration, Instant};

use super::stages::{BooleanSystem, Params, SolverCost, SolverVerdict, SystemShape, SystemSolver};
use crate::cryptanalysis::pq_groebner_f2::{groebner_basis_f2_within, solve_system_f2, F2BoolPoly};

/// The largest system `exhaustive` and the model-enumerating solvers
/// will attempt, since they are `2^n` in the variable count.
const ENUMERATION_CAP: usize = 26;

fn is_one(gb: &[F2BoolPoly]) -> bool {
    gb.len() == 1 && gb[0].terms.len() == 1 && gb[0].terms[0].degree() == 0
}

/// Solutions of a system by direct evaluation, used both as the
/// `exhaustive` solver and to extract solutions from a basis.
fn points_satisfying(equations: &[F2BoolPoly], n_vars: usize) -> (Vec<u64>, u64) {
    let terms: u64 = equations.iter().map(|p| p.terms.len() as u64).sum();
    let mut out = Vec::new();
    for v in 0..1u64 << n_vars {
        if equations.iter().all(|e| e.eval(v) == 0) {
            out.push(v);
        }
    }
    (out, terms * (1u64 << n_vars))
}

// ── Buchberger over the boolean ring ───────────────────────────────

/// The repository's boolean-ring Buchberger, with the coprime and
/// chain criteria under the normal selection strategy.
pub struct BuchbergerF2;

impl SystemSolver for BuchbergerF2 {
    fn name(&self) -> &str {
        "buchberger-f2"
    }

    fn describe(&self) -> String {
        "Buchberger over F_2[v]/(v²−v), DegRevLex, normal selection, coprime and chain criteria"
            .into()
    }

    fn accepts(&self, shape: &SystemShape) -> bool {
        shape.n_vars <= 64
    }

    fn solve(
        &self,
        system: &BooleanSystem,
        _params: &Params,
        budget: Option<Duration>,
    ) -> (SolverVerdict, SolverCost) {
        let (gb, st) = groebner_basis_f2_within(system.equations.clone(), system.n_vars, budget);
        let mut extra = BTreeMap::new();
        extra.insert("spolys".into(), st.spolys);
        extra.insert("reduction_steps".into(), st.reduction_steps);
        extra.insert("pairs_coprime_skipped".into(), st.pairs_coprime_skipped);
        extra.insert("pairs_chain_skipped".into(), st.pairs_chain_skipped);
        extra.insert("basis_len".into(), st.basis_len);
        extra.insert("pairs_left".into(), st.pairs_left);
        let cost = SolverCost {
            ops: st.mono_ops,
            op_unit: "monomial operations".into(),
            wall_ns: st.wall_ns,
            peak_bytes: st.peak_bytes(),
            degree_reached: Some(st.max_pair_degree),
            solving_degree: Some(st.solving_degree),
            timed_out: st.timed_out,
            extra,
        };
        if st.timed_out {
            return (SolverVerdict::BudgetExceeded, cost);
        }
        if is_one(&gb) {
            return (SolverVerdict::Unsatisfiable, cost);
        }
        if system.n_vars > ENUMERATION_CAP {
            // The basis is correct but this module extracts solutions by
            // enumeration, which is not affordable here.  Saying so is
            // better than reporting no solutions.
            return (SolverVerdict::BudgetExceeded, cost);
        }
        let solutions = solve_system_f2(&gb, system.n_vars);
        if solutions.is_empty() {
            (SolverVerdict::Unsatisfiable, cost)
        } else {
            (SolverVerdict::Solved(solutions), cost)
        }
    }
}

// ── XL ─────────────────────────────────────────────────────────────

/// XL (Courtois–Klimov–Patarin–Shamir): multiply the system out to a
/// degree, then linearise.  Scales better than Buchberger on the wider
/// boolean systems and worse on the narrow ones, which is the whole
/// reason to be able to switch between them.
pub struct XlF2;

impl SystemSolver for XlF2 {
    fn name(&self) -> &str {
        "xl-f2"
    }

    fn describe(&self) -> String {
        "XL over F_2: multiply to a degree bound, linearise, solve".into()
    }

    fn accepts(&self, shape: &SystemShape) -> bool {
        shape.n_vars <= ENUMERATION_CAP
    }

    fn solve(
        &self,
        system: &BooleanSystem,
        _params: &Params,
        _budget: Option<Duration>,
    ) -> (SolverVerdict, SolverCost) {
        // The XL implementation has no interrupt hook, so the budget is
        // enforced by `accepts` and by the caller's cell sizing rather
        // than mid-run.  That is stated rather than hidden: a long XL
        // call will run to completion.
        let started = Instant::now();
        let solutions =
            crate::cryptanalysis::pq_xl::boolean_xl_solve(system.equations.clone(), system.n_vars);
        let terms: u64 = system.equations.iter().map(|p| p.terms.len() as u64).sum();
        let cost = SolverCost {
            // XL's work is dominated by the linearised elimination; the
            // term count times the expansion is the honest proxy this
            // implementation can report exactly.
            ops: terms * (system.n_vars as u64).pow(2),
            op_unit: "monomial operations (modelled)".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: None,
            solving_degree: None,
            timed_out: false,
            extra: BTreeMap::new(),
        };
        if solutions.is_empty() {
            (SolverVerdict::Unsatisfiable, cost)
        } else {
            (SolverVerdict::Solved(solutions), cost)
        }
    }
}

// ── CDCL ───────────────────────────────────────────────────────────

/// The repository's CDCL solver, with one Tseitin auxiliary per
/// monomial of degree at least two and the equations as parity rows.
pub struct SatCdcl;

impl SystemSolver for SatCdcl {
    fn name(&self) -> &str {
        "sat-cdcl"
    }

    fn describe(&self) -> String {
        "CDCL with Tseitin monomial definitions and native parity rows".into()
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "sat_conflict_budget",
            "conflicts before the solver gives up (default 200000)",
        )]
    }

    fn accepts(&self, shape: &SystemShape) -> bool {
        shape.n_vars <= 64
    }

    fn solve(
        &self,
        system: &BooleanSystem,
        params: &Params,
        _budget: Option<Duration>,
    ) -> (SolverVerdict, SolverCost) {
        use crate::cryptanalysis::sat::SolveResult;
        let started = Instant::now();
        let mut enc = crate::cryptanalysis::semaev_sat::encode_boolean_system(
            system.n_vars,
            &system.equations,
            &[],
        );
        enc.solver.conflict_budget = params
            .u64_or("sat_conflict_budget", 200_000)
            .unwrap_or(200_000);
        let verdict = if enc.trivially_unsat {
            SolveResult::Unsat
        } else {
            enc.solver.solve()
        };
        let st = &enc.solver.stats;
        let mut extra = BTreeMap::new();
        extra.insert("decisions".into(), st.decisions);
        extra.insert("propagations".into(), st.propagations);
        extra.insert("restarts".into(), st.restarts);
        extra.insert("learnt_clauses".into(), st.learnt_clauses);
        extra.insert("xor_propagations".into(), st.xor_propagations);
        let cost = SolverCost {
            ops: st.conflicts,
            op_unit: "conflicts".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: None,
            solving_degree: None,
            timed_out: matches!(verdict, SolveResult::Unknown),
            extra,
        };
        match verdict {
            SolveResult::Sat => (SolverVerdict::Solved(vec![enc.model_assignment()]), cost),
            SolveResult::Unsat => (SolverVerdict::Unsatisfiable, cost),
            SolveResult::Unknown => (SolverVerdict::BudgetExceeded, cost),
        }
    }
}

// ── Exhaustive search, the reference ───────────────────────────────

/// Evaluate every equation at every point of `{0,1}^n`.
///
/// This is the **reference** every other engine is measured against:
/// the best algorithm that already solves the problem, in the same
/// unit, on the same instance.  Its cost is exactly
/// `2^n · Σ_i |terms_i|` and needs no calibration, which is what makes
/// it a usable denominator.
pub struct Exhaustive;

impl SystemSolver for Exhaustive {
    fn name(&self) -> &str {
        "exhaustive"
    }

    fn describe(&self) -> String {
        "evaluate every equation at every point of {0,1}^n — the reference, not a strawman".into()
    }

    fn accepts(&self, shape: &SystemShape) -> bool {
        shape.n_vars <= ENUMERATION_CAP
    }

    fn solve(
        &self,
        system: &BooleanSystem,
        _params: &Params,
        _budget: Option<Duration>,
    ) -> (SolverVerdict, SolverCost) {
        let started = Instant::now();
        let (solutions, ops) = points_satisfying(&system.equations, system.n_vars);
        let cost = SolverCost {
            ops,
            op_unit: "monomial tests".into(),
            wall_ns: started.elapsed().as_nanos() as u64,
            peak_bytes: 0,
            degree_reached: None,
            solving_degree: None,
            timed_out: false,
            extra: BTreeMap::new(),
        };
        if solutions.is_empty() {
            (SolverVerdict::Unsatisfiable, cost)
        } else {
            (SolverVerdict::Solved(solutions), cost)
        }
    }
}

/// Every solver the framework knows, by name.
pub fn solver_registry() -> Vec<Box<dyn SystemSolver>> {
    vec![
        Box::new(BuchbergerF2),
        Box::new(XlF2),
        Box::new(SatCdcl),
        Box::new(Exhaustive),
    ]
}

/// Look one up.
pub fn solver_by_name(name: &str) -> Result<Box<dyn SystemSolver>, String> {
    solver_registry()
        .into_iter()
        .find(|s| s.name() == name)
        .ok_or_else(|| {
            let known: Vec<&str> = vec!["buchberger-f2", "xl-f2", "sat-cdcl", "exhaustive"];
            format!("unknown solver `{name}`; known: {}", known.join(", "))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

    /// `v0 + v1`, `v0·v1 + v1` — solutions are `v0 = v1 = 0` and
    /// `v0 = v1 = 1`, i.e. masks 0 and 3.
    fn fixture() -> BooleanSystem {
        let m = F2BoolMono::from_mask;
        BooleanSystem {
            equations: vec![
                F2BoolPoly::from_monos(vec![m(1), m(2)], 2),
                F2BoolPoly::from_monos(vec![m(3), m(2)], 2),
            ],
            n_vars: 2,
        }
    }

    /// **Every engine must agree on the answer.**  A framework whose
    /// plug points disagree is not measuring the same problem, and the
    /// column that compares them would be meaningless.
    #[test]
    fn every_solver_finds_the_same_solution_set() {
        let sys = fixture();
        let reference: std::collections::BTreeSet<u64> =
            match Exhaustive.solve(&sys, &Params::default(), None).0 {
                SolverVerdict::Solved(v) => v.into_iter().collect(),
                other => panic!("the reference must solve the fixture, got {other:?}"),
            };
        assert_eq!(reference, [0u64, 3].into_iter().collect());

        for solver in solver_registry() {
            // The CDCL solver returns one model rather than all of
            // them, so it is checked for membership, not equality.
            let (verdict, cost) = solver.solve(&sys, &Params::default(), None);
            match verdict {
                SolverVerdict::Solved(found) => {
                    assert!(!found.is_empty(), "{}: solved with no solutions", solver.name());
                    for f in &found {
                        assert!(
                            reference.contains(f),
                            "{} returned {f}, which is not a solution",
                            solver.name()
                        );
                    }
                    if solver.name() != "sat-cdcl" {
                        let got: std::collections::BTreeSet<u64> = found.into_iter().collect();
                        assert_eq!(got, reference, "{} missed a solution", solver.name());
                    }
                }
                other => panic!("{} said {other:?} on a satisfiable system", solver.name()),
            }
            assert!(
                !cost.op_unit.is_empty(),
                "{} must name the unit it counts",
                solver.name()
            );
        }
    }

    /// An unsatisfiable system must be reported as unsatisfiable by
    /// every engine, never as a budget failure and never as solved.
    #[test]
    fn every_solver_refutes_an_unsatisfiable_system() {
        let m = F2BoolMono::from_mask;
        // v0 = 0 and v0 = 1 at once.
        let sys = BooleanSystem {
            equations: vec![
                F2BoolPoly::from_monos(vec![m(1)], 1),
                F2BoolPoly::from_monos(vec![m(1), m(0)], 1),
            ],
            n_vars: 1,
        };
        for solver in solver_registry() {
            let (verdict, _) = solver.solve(&sys, &Params::default(), None);
            assert!(
                matches!(verdict, SolverVerdict::Unsatisfiable),
                "{} said {verdict:?} on an unsatisfiable system",
                solver.name()
            );
        }
    }

    #[test]
    fn the_registry_resolves_every_name_it_advertises() {
        for solver in solver_registry() {
            let name = solver.name().to_string();
            assert_eq!(solver_by_name(&name).unwrap().name(), name);
        }
        assert!(solver_by_name("f5").is_err());
    }
}
