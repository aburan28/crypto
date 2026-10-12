//! The CDCL shell: CNF layout, static clauses, the algebraic theory (the leaf
//! as a lazy-clause propagator that fires once two of the three summands are
//! assigned) and the arms that drive it.
//!
//! **Variables** (1-based): the bits of `X₁, X₂, X₃`, summand-major, most
//! significant bit first (`1 .. 3ℓ`); one conflict variable `d` (`3ℓ + 1`);
//! the prefix-equal variables `E₁ … E_{ℓ−1}` of the symmetry-breaking
//! constraint (`3ℓ + 2 …`).
//!
//! **Static clauses** (original, not learned): `X₁ ≤ X₂` (the same sorted pair
//! space as the plain loop), and for a candidate set `S ≠ V` one clause per
//! value outside `S` and summand excluding it (the plain loop's candidate set
//! is the rational base points; the shell is given the same set).
//!
//! **Theory.**  When both bit-vectors `X₁, X₂` are assigned the leaf is asked
//! for the third abscissae.  No third in `V`: the explanation is the guard
//! `g` (the negation of the 2ℓ assigned pair literals), installed as
//! `g ∨ d` and `g ∨ ¬d`: the solver propagates `d`, meets the second clause
//! falsified, and analyses an ordinary conflict whose learned clause is the
//! guard or something shorter.  (Installing the bare guard takes the library's
//! conservative path, which backjumps to level 0 and counts no conflict; that
//! encoding is available as `--encoding plain`.)  Some third `t`: the clauses
//! `g ∨ (X₃ = t)` (per bit; all combinations of bit positions for several
//! thirds).

use crate::acct::PairAcct;
use crate::dynorder::{DynStats, MvcDecider, XGraph};
use crate::leaf::Leaf;
use crate::sat_hooked::{self, ShellHooks};
use crypto_lib::cryptanalysis::sat::{self, Lit};
use std::collections::HashSet;

#[derive(Clone, Copy, Debug)]
pub struct Layout {
    pub l: u32,
}

impl Layout {
    pub fn new(l: u32) -> Self {
        Self { l }
    }
    /// SAT variable of bit `j` (0 = least significant) of summand `s`.
    pub fn xv(&self, s: usize, j: u32) -> u32 {
        1 + s as u32 * self.l + (self.l - 1 - j)
    }
    pub fn d(&self) -> u32 {
        3 * self.l + 1
    }
    /// `E_k`, `1 ≤ k ≤ ℓ − 1`.
    pub fn e(&self, k: u32) -> u32 {
        3 * self.l + 1 + k
    }
    pub fn n_vars(&self) -> u32 {
        3 * self.l + self.l.max(1)
    }
    /// The pair variables `X₁, X₂`.
    pub fn pair_vars(&self) -> Vec<u32> {
        (1..=2 * self.l).collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    /// `g ∨ d`, `g ∨ ¬d`: an ordinary conflict, analysed.
    Aux,
    /// The bare guard `g`: the library's level-0 path.
    Plain,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    /// The library solver as shipped: priority on the pair bits, VSIDS within, phase 0 first.
    Native,
    /// The same solver code with `decide` naming the lexicographically first unassigned bit, 0 first.
    Lex,
    /// The same with the re-evaluated minimum vertex cover rule of `dynorder`.
    Mvc,
}

/// Which solver code runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Engine {
    Library,
    Hooked,
}

/// Solver counters, the same for both solver copies.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Stats {
    pub decisions: u64,
    pub conflicts: u64,
    pub restarts: u64,
    pub propagations: u64,
    pub learnt_clauses: u64,
    pub clause_visits: u64,
    pub clause_lit_visits: u64,
    pub analyze_lit_visits: u64,
    pub learnt_lits_raw: u64,
    pub learnt_lits_kept: u64,
    pub max_level: u64,
    pub conflict_level_sum: u64,
}

impl From<sat::SolverStats> for Stats {
    fn from(s: sat::SolverStats) -> Self {
        Stats {
            decisions: s.decisions,
            conflicts: s.conflicts,
            restarts: s.restarts,
            propagations: s.propagations,
            learnt_clauses: s.learnt_clauses,
            clause_visits: s.clause_visits,
            clause_lit_visits: s.clause_lit_visits,
            analyze_lit_visits: s.analyze_lit_visits,
            learnt_lits_raw: s.learnt_lits_raw,
            learnt_lits_kept: s.learnt_lits_kept,
            max_level: s.max_level,
            conflict_level_sum: s.conflict_level_sum,
        }
    }
}

impl From<sat_hooked::SolverStats> for Stats {
    fn from(s: sat_hooked::SolverStats) -> Self {
        Stats {
            decisions: s.decisions,
            conflicts: s.conflicts,
            restarts: s.restarts,
            propagations: s.propagations,
            learnt_clauses: s.learnt_clauses,
            clause_visits: s.clause_visits,
            clause_lit_visits: s.clause_lit_visits,
            analyze_lit_visits: s.analyze_lit_visits,
            learnt_lits_raw: s.learnt_lits_raw,
            learnt_lits_kept: s.learnt_lits_kept,
            max_level: s.max_level,
            conflict_level_sum: s.conflict_level_sum,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    Sat,
    Unsat,
    Unknown,
}

/// The two solver copies behind one interface.
pub trait Sink {
    fn add(&mut self, c: Vec<Lit>) -> bool;
}
impl Sink for sat::Solver {
    fn add(&mut self, c: Vec<Lit>) -> bool {
        self.add_clause(c)
    }
}
impl Sink for sat_hooked::Solver {
    fn add(&mut self, c: Vec<Lit>) -> bool {
        self.add_clause(c)
    }
}

/// Static clauses of the shell for the candidate set `space`.
pub fn static_clauses(lay: Layout, space: &[u64]) -> Vec<Vec<Lit>> {
    let l = lay.l;
    let mut out: Vec<Vec<Lit>> = Vec::new();
    // X₁ ≤ X₂, most significant bit first; E_k: the first k bits are equal.
    for k in 1..=l {
        let j = l - k;
        let (a, b) = (lay.xv(0, j) as Lit, lay.xv(1, j) as Lit);
        let mut leq = vec![-a, b];
        if k > 1 {
            leq.push(-(lay.e(k - 1) as Lit));
        }
        out.push(leq);
        if k < l {
            let e = lay.e(k) as Lit;
            if k > 1 {
                out.push(vec![-e, lay.e(k - 1) as Lit]);
            }
            out.push(vec![-e, -a, b]);
            out.push(vec![-e, a, -b]);
            let mut c1 = vec![e, a, b];
            let mut c2 = vec![e, -a, -b];
            if k > 1 {
                c1.push(-(lay.e(k - 1) as Lit));
                c2.push(-(lay.e(k - 1) as Lit));
            }
            out.push(c1);
            out.push(c2);
        }
    }
    // Values outside the candidate set, per summand.
    if space.len() < (1usize << l) {
        let mut inside = vec![false; 1usize << l];
        for &x in space {
            inside[x as usize] = true;
        }
        for v in 0..(1u64 << l) {
            if inside[v as usize] {
                continue;
            }
            for s in 0..3 {
                out.push(
                    (0..l)
                        .map(|j| {
                            let var = lay.xv(s, j) as Lit;
                            if v >> j & 1 == 1 {
                                -var
                            } else {
                                var
                            }
                        })
                        .collect(),
                );
            }
        }
    }
    out
}

fn summand_value(lay: Layout, asg: &[Option<bool>], s: usize) -> Option<u64> {
    let mut v = 0u64;
    for j in 0..lay.l {
        match asg[(lay.xv(s, j) - 1) as usize] {
            Some(true) => v |= 1 << j,
            Some(false) => {}
            None => return None,
        }
    }
    Some(v)
}

/// The theory and its counters.
pub struct ShellState<L: Leaf> {
    pub leaf: L,
    pub lay: Layout,
    pub encoding: Encoding,
    pub acct: PairAcct,
    installed: HashSet<(u64, u64)>,
    cands: Vec<u64>,
    /// Fixpoints at which both pair summands were assigned.
    pub frontier_calls: u64,
    pub leaf_visits: u64,
    pub refuted: u64,
    pub feasible: u64,
    pub clauses_returned: u64,
    pub literals_returned: u64,
    /// Leaf visits when the first feasible pair was found.
    pub first_feasible_at: Option<u64>,
    /// Scans of the assignment by the decider (words examined).
    pub decide_scans: u64,
}

impl<L: Leaf> ShellState<L> {
    pub fn new(leaf: L, lay: Layout, encoding: Encoding) -> Self {
        let acct = PairAcct::new(lay, leaf.space());
        Self {
            leaf,
            lay,
            encoding,
            acct,
            installed: HashSet::new(),
            cands: Vec::new(),
            frontier_calls: 0,
            leaf_visits: 0,
            refuted: 0,
            feasible: 0,
            clauses_returned: 0,
            literals_returned: 0,
            first_feasible_at: None,
            decide_scans: 0,
        }
    }

    pub fn on_fixpoint(&mut self, asg: &[Option<bool>]) -> Option<Vec<Vec<Lit>>> {
        let lay = self.lay;
        let a = summand_value(lay, asg, 0)?;
        let b = summand_value(lay, asg, 1)?;
        self.frontier_calls += 1;
        if self.installed.contains(&(a, b)) {
            return None;
        }
        debug_assert!(a <= b, "symmetry breaking admits a > b");
        let mut cands = std::mem::take(&mut self.cands);
        self.leaf.eval(a, b, &mut cands);
        self.leaf_visits += 1;
        self.acct.on_eval(a, b);
        let guard: Vec<Lit> = (0..2usize)
            .flat_map(|s| {
                (0..lay.l).map(move |j| {
                    let var = lay.xv(s, j) as Lit;
                    let bit = (if s == 0 { a } else { b }) >> j & 1 == 1;
                    if bit {
                        -var
                    } else {
                        var
                    }
                })
            })
            .collect();
        let clauses: Vec<Vec<Lit>> = if cands.is_empty() {
            self.refuted += 1;
            match self.encoding {
                Encoding::Aux => {
                    let d = lay.d() as Lit;
                    let mut c1 = guard.clone();
                    c1.push(d);
                    let mut c2 = guard;
                    c2.push(-d);
                    vec![c1, c2]
                }
                Encoding::Plain => {
                    self.acct.on_theory_exclusion(a, b);
                    vec![guard]
                }
            }
        } else {
            self.feasible += 1;
            self.first_feasible_at.get_or_insert(self.leaf_visits);
            self.installed.insert((a, b));
            x3_clauses(lay, &guard, &cands)
        };
        self.clauses_returned += clauses.len() as u64;
        self.literals_returned += clauses.iter().map(|c| c.len() as u64).sum::<u64>();
        self.cands = cands;
        Some(clauses)
    }
}

/// `guard ∨ (X₃ ∈ cands)`.
fn x3_clauses(lay: Layout, guard: &[Lit], cands: &[u64]) -> Vec<Vec<Lit>> {
    let l = lay.l;
    let lit = |t: u64, j: u32| -> Lit {
        let var = lay.xv(2, j) as Lit;
        if t >> j & 1 == 1 {
            var
        } else {
            -var
        }
    };
    let mut out = Vec::new();
    let k = cands.len();
    let mut pos = vec![0u32; k];
    loop {
        let mut c = guard.to_vec();
        for (ci, &t) in cands.iter().enumerate() {
            c.push(lit(t, pos[ci]));
        }
        out.push(c);
        let mut i = 0;
        while i < k {
            pos[i] += 1;
            if pos[i] < l {
                break;
            }
            pos[i] = 0;
            i += 1;
        }
        if i == k {
            break;
        }
    }
    out
}

/// Hooks for the example-local solver copy.
pub struct Hooks<'g, L: Leaf> {
    pub st: ShellState<L>,
    pub order: Order,
    pub mvc: Option<MvcDecider<'g>>,
    pub log_learnt: bool,
}

impl<L: Leaf> ShellHooks for Hooks<'_, L> {
    fn theory(&mut self, asg: &[Option<bool>]) -> Option<Vec<Vec<Lit>>> {
        self.st.on_fixpoint(asg)
    }

    fn decide(&mut self, asg: &[Option<bool>]) -> Option<Lit> {
        match self.order {
            Order::Native => None,
            Order::Lex => {
                let nx = 3 * self.st.lay.l as usize;
                for (k, a) in asg.iter().take(nx).enumerate() {
                    self.st.decide_scans += 1;
                    if a.is_none() {
                        return Some(-(k as Lit + 1));
                    }
                }
                None
            }
            Order::Mvc => {
                let pick = self.mvc.as_mut().expect("mvc decider").choose(asg)?;
                Some(-(pick as Lit + 1))
            }
        }
    }

    fn learnt(&mut self, clause: &[Lit]) {
        if self.log_learnt {
            self.st.acct.on_learnt(clause);
        }
    }
}

/// Everything one shell run reports.
pub struct ShellRun {
    pub outcome: Outcome,
    pub witness: Option<[u64; 3]>,
    pub roots: Vec<[u64; 3]>,
    pub stats: Stats,
    pub dyn_stats: Option<DynStats>,
    pub wall_ns: u128,
    pub models: u64,
}

pub struct ShellCfg<'g> {
    pub l: u32,
    pub engine: Engine,
    pub order: Order,
    pub encoding: Encoding,
    pub all_roots: bool,
    pub graph: Option<&'g XGraph>,
    pub log_learnt: bool,
    pub max_conflicts: u64,
}

fn sorted3(a: u64, b: u64, c: u64) -> [u64; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

fn decode(lay: Layout, model: &[bool]) -> [u64; 3] {
    let mut v = [0u64; 3];
    for (s, slot) in v.iter_mut().enumerate() {
        for j in 0..lay.l {
            if model[(lay.xv(s, j) - 1) as usize] {
                *slot |= 1 << j;
            }
        }
    }
    v
}

fn blocking(lay: Layout, model: &[bool]) -> Vec<Lit> {
    (1..=3 * lay.l)
        .map(|v| {
            if model[(v - 1) as usize] {
                -(v as Lit)
            } else {
                v as Lit
            }
        })
        .collect()
}

/// Run the shell on one target.  Returns the run and the final state (for counters).
pub fn run_shell<'g, L: Leaf>(cfg: &ShellCfg<'g>, leaf: L) -> (ShellRun, ShellState<L>) {
    let lay = Layout::new(cfg.l);
    let statics = static_clauses(lay, leaf.space());
    let st = ShellState::new(leaf, lay, cfg.encoding);
    let mvc = match cfg.order {
        Order::Mvc => Some(MvcDecider::new(cfg.graph.expect("--order mvc needs the graph"))),
        _ => None,
    };
    let mut hooks = Hooks {
        st,
        order: cfg.order,
        mvc,
        log_learnt: cfg.log_learnt,
    };
    let pair = lay.pair_vars();
    let t0 = std::time::Instant::now();
    let mut roots = Vec::new();
    let mut witness = None;
    let mut models = 0u64;
    let mut outcome;
    let stats: Stats;
    match cfg.engine {
        Engine::Library => {
            let mut s = sat::Solver::new(lay.n_vars());
            for c in statics {
                s.add(c);
            }
            s.set_branch_priority(&pair);
            s.set_all_saved_phases(false);
            s.conflict_budget = cfg.max_conflicts;
            loop {
                let mut theory = |asg: &[Option<bool>]| hooks.st.on_fixpoint(asg);
                let r = s.solve_with_lazy_clauses(&pair, &mut theory);
                outcome = match r {
                    sat::SolveResult::Sat => Outcome::Sat,
                    sat::SolveResult::Unsat => Outcome::Unsat,
                    sat::SolveResult::Unknown => Outcome::Unknown,
                };
                if outcome != Outcome::Sat {
                    break;
                }
                let model = s.model();
                let [a, b, c] = decode(lay, &model);
                models += 1;
                let tr = sorted3(a, b, c);
                if witness.is_none() {
                    witness = Some(tr);
                }
                roots.push(tr);
                if !cfg.all_roots {
                    break;
                }
                s.reset_search();
                if !s.add(blocking(lay, &model)) {
                    outcome = Outcome::Unsat;
                    break;
                }
            }
            stats = s.stats.into();
        }
        Engine::Hooked => {
            let mut s = sat_hooked::Solver::new(lay.n_vars());
            for c in statics {
                s.add(c);
            }
            s.set_branch_priority(&pair);
            s.set_all_saved_phases(false);
            s.conflict_budget = cfg.max_conflicts;
            // The library's trigger rule for the native order (the fidelity check); the other
            // orders let the theory test for a complete pair itself, at every fixpoint.
            let trigger: Option<&[u32]> = if cfg.order == Order::Native {
                Some(&pair)
            } else {
                None
            };
            loop {
                let r = s.solve_with_hooks(trigger, &mut hooks);
                outcome = match r {
                    sat_hooked::SolveResult::Sat => Outcome::Sat,
                    sat_hooked::SolveResult::Unsat => Outcome::Unsat,
                    sat_hooked::SolveResult::Unknown => Outcome::Unknown,
                };
                if outcome != Outcome::Sat {
                    break;
                }
                let model = s.model();
                let [a, b, c] = decode(lay, &model);
                models += 1;
                let tr = sorted3(a, b, c);
                if witness.is_none() {
                    witness = Some(tr);
                }
                roots.push(tr);
                if !cfg.all_roots {
                    break;
                }
                s.reset_search();
                if !s.add(blocking(lay, &model)) {
                    outcome = Outcome::Unsat;
                    break;
                }
            }
            stats = s.stats.into();
        }
    }
    let wall_ns = t0.elapsed().as_nanos();
    roots.sort_unstable();
    roots.dedup();
    let dyn_stats = hooks.mvc.as_ref().map(|m| m.stats.clone());
    (
        ShellRun {
            outcome,
            witness,
            roots,
            stats,
            dyn_stats,
            wall_ns,
            models,
        },
        hooks.st,
    )
}

/// The plain loop over a leaf, one pair at a time (A1 for the quartic leaf, A1u for the group leaf).
/// Returns (found, pairs visited, roots).
pub fn plain_loop<L: Leaf>(leaf: &mut L, all: bool) -> (bool, u64, Vec<[u64; 3]>) {
    let space: Vec<u64> = leaf.space().to_vec();
    let mut pairs = 0u64;
    let mut roots = Vec::new();
    let mut found = false;
    let mut out = Vec::new();
    'rows: for (i, &a) in space.iter().enumerate() {
        for &b in &space[i..] {
            pairs += 1;
            leaf.eval(a, b, &mut out);
            for &t in &out {
                found = true;
                roots.push(sorted3(a, b, t));
            }
            if found && !all {
                break 'rows;
            }
        }
    }
    roots.sort_unstable();
    roots.dedup();
    (found, pairs, roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A leaf over a tiny explicit relation, to test the shell against enumeration.
    struct ToyLeaf {
        space: Vec<u64>,
        good: Vec<(u64, u64, u64)>,
        evals: u64,
    }
    impl Leaf for ToyLeaf {
        fn space(&self) -> &[u64] {
            &self.space
        }
        fn eval(&mut self, a: u64, b: u64, out: &mut Vec<u64>) {
            out.clear();
            self.evals += 1;
            for &(x, y, t) in &self.good {
                if (x, y) == (a, b) && !out.contains(&t) {
                    out.push(t);
                }
            }
        }
        fn counts(&self) -> crate::leaf::LeafCounts {
            Default::default()
        }
        fn name(&self) -> &'static str {
            "toy"
        }
    }

    fn toy(l: u32, seed: u64) -> (ToyLeaf, Vec<[u64; 3]>) {
        let mut s = seed | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let span = 1u64 << l;
        let space: Vec<u64> = (0..span).filter(|_| next() % 4 != 0).collect();
        let mut good = Vec::new();
        for (i, &a) in space.iter().enumerate() {
            for &b in &space[i..] {
                if next() % 40 == 0 {
                    good.push((a, b, space[(next() % space.len() as u64) as usize]));
                }
            }
        }
        let mut roots: Vec<[u64; 3]> = good.iter().map(|&(a, b, t)| sorted3(a, b, t)).collect();
        roots.sort_unstable();
        roots.dedup();
        (
            ToyLeaf {
                space,
                good,
                evals: 0,
            },
            roots,
        )
    }

    #[test]
    fn every_arm_finds_exactly_the_roots_of_a_toy_relation() {
        for l in [3u32, 4] {
            for seed in 1..6u64 {
                let (leaf, roots) = toy(l, seed * 77 + l as u64);
                for (engine, order, enc) in [
                    (Engine::Library, Order::Native, Encoding::Aux),
                    (Engine::Library, Order::Native, Encoding::Plain),
                    (Engine::Hooked, Order::Native, Encoding::Aux),
                    (Engine::Hooked, Order::Lex, Encoding::Aux),
                    (Engine::Hooked, Order::Lex, Encoding::Plain),
                ] {
                    let cfg = ShellCfg {
                        l,
                        engine,
                        order,
                        encoding: enc,
                        all_roots: true,
                        graph: None,
                        log_learnt: true,
                        max_conflicts: u64::MAX,
                    };
                    let (run, _st) = run_shell(
                        &cfg,
                        ToyLeaf {
                            space: leaf.space.clone(),
                            good: leaf.good.clone(),
                            evals: 0,
                        },
                    );
                    assert_eq!(run.outcome, Outcome::Unsat, "{engine:?} {order:?} {enc:?}");
                    assert_eq!(run.roots, roots, "{engine:?} {order:?} {enc:?} l={l} seed={seed}");
                }
            }
        }
    }
}
