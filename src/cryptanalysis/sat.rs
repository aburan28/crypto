//! From-scratch CDCL SAT solver for cryptanalytic boolean problems.
//!
//! Modern SAT solving for cryptanalysis goes back to Massacci-Marraro
//! (2000) and exploded after Soos-Nohl-Castelluccia's CryptoMiniSat
//! (2009) added XOR-clause reasoning. The technique is now standard for:
//!
//! - **Algebraic attacks**: solve the polynomial system describing a
//!   reduced-round cipher (see `cryptanalysis::aes::algebraic`).
//! - **Differential trail search**: encode active S-box constraints as
//!   clauses (Mouha-Preneel 2013, Sun et al. 2014).
//! - **Preimage / collision search**: encode reduced hash functions.
//! - **Side-channel template inversion**: encode leakage equations.
//!
//! This module ships a **Conflict-Driven Clause Learning (CDCL)** SAT
//! solver from scratch — no external crate. The implementation
//! follows the standard recipe (MiniSAT/Glucose family):
//!
//! - Two-watched-literal unit propagation,
//! - 1-UIP conflict analysis with clause learning,
//! - Non-chronological backjumping,
//! - VSIDS-like variable activity heuristic with exponential decay,
//! - Phase saving,
//! - Luby-sequence restarts.
//!
//! It is **not competitive** with state-of-the-art solvers like
//! kissat. It's didactic — every step is readable, the data layout
//! matches the textbook, and it solves small cryptanalytic instances
//! (pigeonhole, n-queens, the AES 4-bit S-box system) in milliseconds.
//!
//! ## DIMACS I/O
//!
//! [`parse_dimacs`] parses the standard `.cnf` input format used by
//! every SAT competition since 1992; [`to_dimacs`] emits it. This
//! lets you round-trip with `kissat`, `cadical`, `minisat`, or any
//! external solver while developing.
//!
//! ## Worked example
//!
//! ```
//! use crypto::cryptanalysis::sat::{Solver, SolveResult};
//!
//! let mut s = Solver::new(3);
//! // (x1 ∨ x2) ∧ (¬x1 ∨ x3) ∧ (¬x2 ∨ ¬x3)
//! s.add_clause(vec![1, 2]);
//! s.add_clause(vec![-1, 3]);
//! s.add_clause(vec![-2, -3]);
//! assert!(matches!(s.solve(), SolveResult::Sat));
//! let model = s.model();
//! // Verify the model satisfies all clauses.
//! ```

use std::collections::HashSet;

/// Standard SAT literal encoding: positive integer = positive
/// literal, negative integer = negated literal, `|lit| - 1` = variable
/// index (0-based internally; 1-based in DIMACS).
pub type Lit = i32;

#[inline]
fn var_of(lit: Lit) -> u32 {
    (lit.unsigned_abs() - 1) as u32
}

#[inline]
fn is_neg(lit: Lit) -> bool {
    lit < 0
}

/// Outcome of solving.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveResult {
    Sat,
    Unsat,
    /// Hit conflict / restart budget without resolving.
    Unknown,
}

/// Reason a literal was assigned: either a unit-propagation source
/// clause, or a free decision.
#[derive(Debug, Clone, Copy)]
enum Reason {
    Decision,
    Propagated(usize), // clause index in `clauses`
}

/// CDCL SAT solver.
pub struct Solver {
    /// Number of variables (1-indexed externally).
    n_vars: u32,
    /// Clauses, original + learnt; learnt start at `n_orig_clauses`.
    clauses: Vec<Vec<Lit>>,
    n_orig_clauses: usize,
    /// Per-variable assignment: None = unassigned.
    assignment: Vec<Option<bool>>,
    /// Per-variable decision level when assigned.
    level: Vec<i32>,
    /// Per-variable reason for assignment.
    reason: Vec<Reason>,
    /// Per-variable last-phase, for phase saving.
    saved_phase: Vec<bool>,
    /// Per-variable VSIDS activity.
    activity: Vec<f64>,
    activity_inc: f64,
    activity_decay: f64,
    /// Trail: order in which literals were assigned.
    trail: Vec<Lit>,
    /// Index into `trail` of first unpropagated literal.
    qhead: usize,
    /// `trail_lim[d] = trail.len() when decision level d began`.
    trail_lim: Vec<usize>,
    /// Two-watched-literal scheme: for each literal (encoded as
    /// `2*var + neg` so we can index a `Vec`), a list of clause
    /// indices watching it.
    watches: Vec<Vec<usize>>,
    /// Counter of conflicts since the last restart.
    conflicts_since_restart: u64,
    /// Total conflicts across the whole solve.
    conflicts: u64,
    /// Maximum total conflicts before giving up. `u64::MAX` means
    /// no limit.
    pub conflict_budget: u64,
    /// Set when `add_clause` detects UNSAT (empty clause or
    /// conflict-on-unit).  `solve()` short-circuits to UNSAT when set.
    is_unsat: bool,
    /// Native XOR (parity) constraints, held *outside* the CNF.  Each
    /// row is `(variable-set bitmask, right-hand side)`, meaning
    /// `⊕_{v ∈ mask} x_v = rhs`.  See [`Solver::add_xor`].
    xors: Vec<XorRow>,
    /// Bumped on every assignment and every backjump.  Lets
    /// [`Solver::propagate`] skip the Gauss-Jordan pass when nothing
    /// has moved since the last one.
    epoch: u64,
    /// `epoch` as of the last completed Gauss-Jordan pass.
    xor_epoch: u64,
}

/// One parity constraint: `⊕_{v ∈ mask} x_v = rhs`, with `mask` a
/// bitmask over 0-indexed variables.
#[derive(Clone, Debug)]
struct XorRow {
    mask: Vec<u64>,
    rhs: bool,
}

/// Outcome of one Gauss-Jordan pass over the XOR rows.
enum XorStep {
    /// Nothing new could be derived.
    Fixpoint,
    /// At least one literal was enqueued; re-run clause propagation.
    Propagated,
    /// A row is inconsistent; the payload is the index of a freshly
    /// installed conflict clause, falsified by the current trail.
    Conflict(usize),
}

#[inline]
fn bs_words(n_vars: u32) -> usize {
    (n_vars as usize + 63) / 64
}

#[inline]
fn bs_get(mask: &[u64], v: u32) -> bool {
    (mask[v as usize / 64] >> (v % 64)) & 1 == 1
}

#[inline]
fn bs_flip(mask: &mut [u64], v: u32) {
    mask[v as usize / 64] ^= 1u64 << (v % 64);
}

#[inline]
fn watch_index(lit: Lit) -> usize {
    let v = var_of(lit) as usize;
    2 * v + (if is_neg(lit) { 1 } else { 0 })
}

#[inline]
fn negated(lit: Lit) -> Lit {
    -lit
}

impl Solver {
    /// Create a fresh solver for `n_vars` variables, indexed 1..=n_vars.
    pub fn new(n_vars: u32) -> Self {
        let n = n_vars as usize;
        Solver {
            n_vars,
            clauses: Vec::new(),
            n_orig_clauses: 0,
            assignment: vec![None; n],
            level: vec![-1; n],
            reason: vec![Reason::Decision; n],
            saved_phase: vec![true; n],
            activity: vec![0.0; n],
            activity_inc: 1.0,
            activity_decay: 0.95,
            trail: Vec::new(),
            qhead: 0,
            trail_lim: Vec::new(),
            watches: vec![Vec::new(); 2 * n],
            conflicts_since_restart: 0,
            conflicts: 0,
            conflict_budget: u64::MAX,
            is_unsat: false,
            xors: Vec::new(),
            epoch: 0,
            xor_epoch: u64::MAX,
        }
    }

    /// **Add a native XOR (parity) constraint** `x_{v₁} ⊕ … ⊕ x_{v_k}
    /// = rhs`, with variables given 1-indexed exactly as in
    /// [`Solver::add_clause`].  Returns `false` if the constraint is
    /// trivially unsatisfiable (`0 = 1`).
    ///
    /// XOR rows are held outside the CNF and reasoned about by
    /// Gauss-Jordan elimination in [`Solver::propagate`], rather than
    /// being Tseitin-expanded into clauses.  This is the difference
    /// between polynomial-time and exponential-time handling of a
    /// dense parity constraint: resolution needs exponentially many
    /// steps to refute one (Urquhart 1987), which is why a plain CDCL
    /// stalls on Weil-descended Semaev systems.
    ///
    /// Duplicated variables cancel (`x ⊕ x = 0`), so the caller need
    /// not deduplicate.
    pub fn add_xor(&mut self, vars: &[u32], rhs: bool) -> bool {
        let mut mask = vec![0u64; bs_words(self.n_vars)];
        for &v in vars {
            debug_assert!(v >= 1 && v <= self.n_vars, "xor var {v} out of range");
            bs_flip(&mut mask, v - 1); // duplicates cancel
        }
        if mask.iter().all(|w| *w == 0) {
            // Empty parity: `0 = rhs`.  Satisfiable iff rhs is false.
            if rhs {
                self.is_unsat = true;
                return false;
            }
            return true;
        }
        self.xors.push(XorRow { mask, rhs });
        self.xor_epoch = u64::MAX; // force a pass on the next propagate
        true
    }

    /// Number of native XOR constraints installed.
    pub fn n_xors(&self) -> usize {
        self.xors.len()
    }

    /// Verify a model against the installed XOR rows.  The CNF part is
    /// checked separately by [`check_model`].
    pub fn check_xors(&self, model: &[bool]) -> bool {
        self.xors.iter().all(|row| {
            let mut parity = false;
            for v in 0..self.n_vars {
                if bs_get(&row.mask, v) && model[v as usize] {
                    parity = !parity;
                }
            }
            parity == row.rhs
        })
    }

    /// Add a clause `lits` (DIMACS literal encoding). Returns false if
    /// the empty clause is added (immediately UNSAT).
    pub fn add_clause(&mut self, mut lits: Vec<Lit>) -> bool {
        // Deduplicate and remove tautologies.
        lits.sort_by_key(|&l| (var_of(l), l));
        lits.dedup();
        for w in lits.windows(2) {
            if w[0] == -w[1] {
                return true; // tautology — drop
            }
        }
        match lits.len() {
            0 => {
                self.is_unsat = true;
                false
            }
            1 => {
                // Unit clause: assign now (at level 0).
                let l = lits[0];
                match self.enqueue(l, Reason::Propagated(self.clauses.len())) {
                    Ok(()) => {
                        self.clauses.push(vec![l]);
                        self.n_orig_clauses = self.clauses.len();
                        true
                    }
                    Err(()) => {
                        self.is_unsat = true;
                        false
                    }
                }
            }
            _ => {
                let idx = self.clauses.len();
                self.watches[watch_index(lits[0])].push(idx);
                self.watches[watch_index(lits[1])].push(idx);
                self.clauses.push(lits);
                self.n_orig_clauses = self.clauses.len();
                true
            }
        }
    }

    /// Look up the truth value of a literal under the current trail.
    fn lit_value(&self, lit: Lit) -> Option<bool> {
        let v = var_of(lit) as usize;
        self.assignment[v].map(|b| if is_neg(lit) { !b } else { b })
    }

    /// Assign `lit` to true with the given reason. Returns `Err` if it
    /// conflicts with the current assignment.
    fn enqueue(&mut self, lit: Lit, r: Reason) -> Result<(), ()> {
        match self.lit_value(lit) {
            Some(true) => Ok(()),
            Some(false) => Err(()),
            None => {
                let v = var_of(lit) as usize;
                self.assignment[v] = Some(!is_neg(lit));
                self.level[v] = self.trail_lim.len() as i32;
                self.reason[v] = r;
                self.saved_phase[v] = !is_neg(lit);
                self.trail.push(lit);
                self.epoch += 1;
                Ok(())
            }
        }
    }

    /// Propagate to fixpoint over *both* reasoning engines: watched-
    /// literal unit propagation over the CNF, and Gauss-Jordan
    /// elimination over the native XOR rows.  Each engine can feed the
    /// other, so we alternate until neither derives anything new.
    ///
    /// Returns `Some(clause_idx)` on conflict.  For an XOR conflict the
    /// index points at a freshly installed clause that is falsified by
    /// the current trail, so `analyze()` and `backjump()` handle it with
    /// no special-casing.
    fn propagate(&mut self) -> Option<usize> {
        loop {
            if let Some(c) = self.propagate_clauses() {
                return Some(c);
            }
            if self.xors.is_empty() || self.epoch == self.xor_epoch {
                // Nothing has moved since the last Gauss-Jordan pass.
                return None;
            }
            self.xor_epoch = self.epoch;
            match self.propagate_xors() {
                XorStep::Conflict(idx) => return Some(idx),
                XorStep::Propagated => continue,
                XorStep::Fixpoint => return None,
            }
        }
    }

    /// One Gauss-Jordan pass over the XOR rows under the current trail.
    ///
    /// Each row carries the *full* variable set of its combination, and
    /// its right-hand side is the combination's original rhs; the value
    /// implied by the current assignment is therefore
    /// `rhs ⊕ parity(mask ∩ assigned-true)`.  Because both mask and rhs
    /// XOR when two rows are combined, that relation survives
    /// elimination — which is what lets us read a reason clause straight
    /// off a reduced row.
    fn propagate_xors(&mut self) -> XorStep {
        let words = bs_words(self.n_vars);
        let mut rows = self.xors.clone();

        // Forward elimination.  Processing rows in order and clearing
        // each chosen pivot from *every* other row keeps the invariant
        // that a pivot lives in exactly one row, so a single pass
        // suffices.  Any pivot already claimed by an earlier row has
        // been eliminated from this one by the time we reach it.
        for i in 0..rows.len() {
            let pivot = match self.lowest_unassigned(&rows[i].mask) {
                Some(p) => p,
                None => continue, // fully assigned: checked below
            };
            let (src_rhs, src_mask) = (rows[i].rhs, rows[i].mask.clone());
            for j in 0..rows.len() {
                if j == i || !bs_get(&rows[j].mask, pivot) {
                    continue;
                }
                rows[j].rhs ^= src_rhs;
                for w in 0..words {
                    rows[j].mask[w] ^= src_mask[w];
                }
            }
        }

        // Read off conflicts and unit implications.
        let mut propagated = false;
        for row in &rows {
            let mut unassigned: Option<u32> = None;
            let mut extra_unassigned = false;
            let mut parity = row.rhs;
            for v in 0..self.n_vars {
                if !bs_get(&row.mask, v) {
                    continue;
                }
                match self.assignment[v as usize] {
                    Some(true) => parity ^= true,
                    Some(false) => {}
                    None => {
                        if unassigned.is_some() {
                            extra_unassigned = true;
                            break;
                        }
                        unassigned = Some(v);
                    }
                }
            }
            if extra_unassigned {
                continue; // under-determined; nothing to derive yet
            }
            match unassigned {
                None => {
                    // Fully assigned.  `parity` is the residual: it must
                    // be 0, or the row is violated.
                    if parity {
                        let clause = self.xor_reason_clause(&row.mask, None);
                        let idx = self.clauses.len();
                        self.clauses.push(clause);
                        return XorStep::Conflict(idx);
                    }
                }
                Some(x) => {
                    // The row forces `x = parity`.
                    let lit = if parity { (x + 1) as Lit } else { -((x + 1) as Lit) };
                    if self.lit_value(lit) == Some(true) {
                        continue; // already implied
                    }
                    let clause = self.xor_reason_clause(&row.mask, Some(lit));
                    let idx = self.clauses.len();
                    self.clauses.push(clause);
                    if self.enqueue(lit, Reason::Propagated(idx)).is_err() {
                        return XorStep::Conflict(idx);
                    }
                    propagated = true;
                }
            }
        }

        if propagated {
            XorStep::Propagated
        } else {
            XorStep::Fixpoint
        }
    }

    /// Lowest-numbered unassigned variable in `mask`, if any.
    fn lowest_unassigned(&self, mask: &[u64]) -> Option<u32> {
        for (w, &word) in mask.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let b = bits.trailing_zeros();
                bits &= bits - 1;
                let v = (w * 64) as u32 + b;
                if v < self.n_vars && self.assignment[v as usize].is_none() {
                    return Some(v);
                }
            }
        }
        None
    }

    /// Build the clause witnessing what a reduced XOR row implies:
    /// every assigned variable of the row appears negated-as-assigned,
    /// so the clause is false under the current trail except for
    /// `implied` (absent for a conflict clause, which is wholly false).
    fn xor_reason_clause(&self, mask: &[u64], implied: Option<Lit>) -> Vec<Lit> {
        let mut clause = Vec::new();
        if let Some(l) = implied {
            clause.push(l);
        }
        let implied_var = implied.map(var_of);
        for v in 0..self.n_vars {
            if !bs_get(mask, v) || Some(v) == implied_var {
                continue;
            }
            match self.assignment[v as usize] {
                Some(true) => clause.push(-((v + 1) as Lit)),
                Some(false) => clause.push((v + 1) as Lit),
                None => debug_assert!(false, "reason clause over an unassigned variable"),
            }
        }
        clause
    }

    /// Unit propagation over the CNF. Returns `Some(clause_idx)` on conflict.
    fn propagate_clauses(&mut self) -> Option<usize> {
        while self.qhead < self.trail.len() {
            let lit = self.trail[self.qhead];
            self.qhead += 1;
            // Iterate clauses watching `¬lit` (the watcher becomes
            // false when `lit` is set true).
            let wi = watch_index(-lit);
            let mut i = 0usize;
            // Move all watchers out so we can mutate `self`.
            let watchers = std::mem::take(&mut self.watches[wi]);
            let mut new_watchers: Vec<usize> = Vec::with_capacity(watchers.len());
            'outer: while i < watchers.len() {
                let cidx = watchers[i];
                i += 1;
                // Ensure clause[0] is the false watcher (the one we're
                // looking to replace).  The watch entry tells us the
                // false watcher is one of clause[0]/clause[1]; swap so
                // it's at index 0.
                {
                    let clause = &mut self.clauses[cidx];
                    if clause[0] != -lit {
                        clause.swap(0, 1);
                    }
                }
                // Read the "other" watcher value (no clause borrow held).
                let other = self.clauses[cidx][1];
                let other_value = self.lit_value(other);
                if other_value == Some(true) {
                    new_watchers.push(cidx);
                    continue;
                }
                // Find a new watch — scan positions 2..len under
                // narrow re-borrow per iteration.
                let clen = self.clauses[cidx].len();
                let mut found_new_watch = false;
                for k in 2..clen {
                    let l = self.clauses[cidx][k];
                    if self.lit_value(l) != Some(false) {
                        // Use `l` as the new watch.
                        self.clauses[cidx].swap(0, k);
                        self.watches[watch_index(l)].push(cidx);
                        found_new_watch = true;
                        break;
                    }
                }
                if found_new_watch {
                    continue 'outer;
                }
                // No new watch: clause is unit or conflict. Restore
                // the canonical layout (false-watcher at [0]).
                self.clauses[cidx][0] = -lit;
                new_watchers.push(cidx);
                match other_value {
                    None => {
                        // Unit propagate `other`.
                        if self.enqueue(other, Reason::Propagated(cidx)).is_err() {
                            // Shouldn't happen given other was None.
                            self.watches[wi].extend(new_watchers);
                            self.watches[wi].extend_from_slice(&watchers[i..]);
                            return Some(cidx);
                        }
                    }
                    Some(false) => {
                        // Conflict: all literals are false.
                        self.watches[wi].extend(new_watchers);
                        self.watches[wi].extend_from_slice(&watchers[i..]);
                        return Some(cidx);
                    }
                    Some(true) => unreachable!(),
                }
            }
            self.watches[wi] = new_watchers;
        }
        None
    }

    /// 1-UIP conflict analysis. Returns `(learnt_clause, backjump_level)`.
    fn analyze(&mut self, conflict_idx: usize) -> (Vec<Lit>, i32) {
        let current_level = self.trail_lim.len() as i32;
        let mut learnt: Vec<Lit> = Vec::new();
        let mut seen: Vec<bool> = vec![false; self.n_vars as usize];
        let mut counter = 0i32;
        let mut p: Lit = 0;
        let mut p_reason_idx = conflict_idx;
        let mut trail_pos = self.trail.len();

        loop {
            // Collect literals from current reason clause.
            let clause: Vec<Lit> = self.clauses[p_reason_idx].clone();
            for &q in &clause {
                if p != 0 && q == p {
                    continue;
                }
                let v = var_of(q) as usize;
                if !seen[v] && self.level[v] >= 0 {
                    seen[v] = true;
                    // Bump activity.
                    self.activity[v] += self.activity_inc;
                    if self.activity[v] > 1e100 {
                        for a in self.activity.iter_mut() {
                            *a *= 1e-100;
                        }
                        self.activity_inc *= 1e-100;
                    }
                    if self.level[v] >= current_level {
                        counter += 1;
                    } else {
                        learnt.push(q);
                    }
                }
            }
            // Find the next literal to resolve on — walk back the trail.
            while trail_pos > 0 {
                trail_pos -= 1;
                let l = self.trail[trail_pos];
                if seen[var_of(l) as usize] {
                    p = l;
                    break;
                }
            }
            let v = var_of(p) as usize;
            seen[v] = false;
            counter -= 1;
            if counter <= 0 {
                break;
            }
            // The reason of p must be a propagation (not a decision).
            match self.reason[v] {
                Reason::Propagated(idx) => p_reason_idx = idx,
                Reason::Decision => break,
            }
        }
        // Asserting literal: ¬p.
        learnt.insert(0, -p);
        // Backjump level: the second-highest level in the clause.
        let mut bj_level = 0;
        if learnt.len() > 1 {
            // Find max level among lits[1..].
            let mut max_i = 1;
            for i in 2..learnt.len() {
                if self.level[var_of(learnt[i]) as usize]
                    > self.level[var_of(learnt[max_i]) as usize]
                {
                    max_i = i;
                }
            }
            learnt.swap(1, max_i);
            bj_level = self.level[var_of(learnt[1]) as usize];
        }
        // Decay activity.
        self.activity_inc /= self.activity_decay;
        (learnt, bj_level)
    }

    /// Undo all assignments above `level`.
    fn backjump(&mut self, level: i32) {
        let level = level.max(0) as usize;
        if self.trail_lim.len() <= level {
            return;
        }
        let target = self.trail_lim[level];
        while self.trail.len() > target {
            let l = self.trail.pop().unwrap();
            let v = var_of(l) as usize;
            self.assignment[v] = None;
            self.level[v] = -1;
        }
        self.trail_lim.truncate(level);
        self.qhead = target;
        self.epoch += 1;
    }

    /// Pick an unassigned variable with the highest activity. Returns
    /// `None` if all variables are assigned.
    fn pick_branching_variable(&self) -> Option<u32> {
        let mut best: Option<(u32, f64)> = None;
        for v in 0..self.n_vars as usize {
            if self.assignment[v].is_some() {
                continue;
            }
            let a = self.activity[v];
            if best.map_or(true, |(_, ba)| a > ba) {
                best = Some((v as u32, a));
            }
        }
        best.map(|(v, _)| v)
    }

    /// Main solve loop. Runs until SAT/UNSAT or conflict budget hits.
    pub fn solve(&mut self) -> SolveResult {
        if self.is_unsat {
            return SolveResult::Unsat;
        }
        // Initial propagation at level 0 catches trivial UNSAT from
        // unit clauses.
        if self.propagate().is_some() {
            return SolveResult::Unsat;
        }
        let mut luby_index = 1u64;
        let mut restart_limit = 100u64 * luby(luby_index);
        loop {
            if let Some(conflict_idx) = self.propagate() {
                self.conflicts += 1;
                self.conflicts_since_restart += 1;
                if self.conflicts >= self.conflict_budget {
                    return SolveResult::Unknown;
                }
                if self.trail_lim.is_empty() {
                    return SolveResult::Unsat;
                }
                let (learnt, bj_level) = self.analyze(conflict_idx);
                self.backjump(bj_level);
                // Install the learnt clause.
                if learnt.len() == 1 {
                    // It's a unit at the backjump level (which is 0).
                    let _ = self.enqueue(learnt[0], Reason::Propagated(self.clauses.len()));
                    self.clauses.push(learnt);
                } else {
                    let idx = self.clauses.len();
                    let l0 = learnt[0];
                    let l1 = learnt[1];
                    self.clauses.push(learnt);
                    self.watches[watch_index(l0)].push(idx);
                    self.watches[watch_index(l1)].push(idx);
                    // The first lit is the asserting one — it's
                    // implied at the new (lower) level.
                    let _ = self.enqueue(l0, Reason::Propagated(idx));
                }
                // Restart?
                if self.conflicts_since_restart >= restart_limit {
                    self.backjump(0);
                    self.conflicts_since_restart = 0;
                    luby_index += 1;
                    restart_limit = 100u64 * luby(luby_index);
                }
            } else {
                // No conflict — pick a new variable.
                match self.pick_branching_variable() {
                    None => return SolveResult::Sat,
                    Some(v) => {
                        self.trail_lim.push(self.trail.len());
                        let lit = if self.saved_phase[v as usize] {
                            (v + 1) as i32
                        } else {
                            -((v + 1) as i32)
                        };
                        let _ = self.enqueue(lit, Reason::Decision);
                    }
                }
            }
        }
    }

    /// Get the satisfying assignment after a successful `solve()`.
    pub fn model(&self) -> Vec<bool> {
        self.assignment.iter().map(|a| a.unwrap_or(false)).collect()
    }

    /// Number of variables.
    pub fn n_vars(&self) -> u32 {
        self.n_vars
    }

    /// Total number of (original + learnt) clauses.
    /// Conflicts encountered by the most recent [`Self::solve`].
    ///
    /// A machine-independent measure of search effort: unlike wall
    /// clock it is comparable across runs and machines, which is what
    /// a benchmark wants when comparing encodings.
    pub fn conflicts(&self) -> u64 {
        self.conflicts
    }

    pub fn n_clauses(&self) -> usize {
        self.clauses.len()
    }
}

/// Luby sequence `1, 1, 2, 1, 1, 2, 4, 1, 1, 2, 1, 1, 2, 4, 8, …`
/// (Luby-Sinclair-Zuckerman 1993). Used for restart timing.
/// Iterative — the natural recursion grows logarithmically but on
/// pathological inputs we don't want to risk the stack.
fn luby(i: u64) -> u64 {
    // i is 1-indexed in this codebase; the sequence starts at i = 1.
    // Luby recurrence (Luby-Sinclair-Zuckerman 1993):
    //   t_i = 2^{k-1}            if i = 2^k - 1
    //        t_{i - 2^{k-1} + 1}  otherwise (with k = ⌊log₂ i⌋ + 1)
    let mut n = i;
    loop {
        if n == 0 {
            return 1; // guard; shouldn't happen with i >= 1
        }
        // Find the smallest k with 2^k > n  (i.e. 2^k >= n + 1).
        let mut k: u64 = 0;
        while (1u64 << k) < n + 1 {
            k += 1;
        }
        // Now (1 << k) >= n + 1.  Two cases:
        if n + 1 == 1u64 << k {
            // n == 2^k - 1: base case, return 2^{k-1} (= 1 when k = 0).
            return if k == 0 { 1 } else { 1u64 << (k - 1) };
        }
        // Otherwise, recurse on n - (2^{k-1} - 1).
        n -= (1u64 << (k - 1)) - 1;
    }
}

/// Parse a DIMACS CNF string.
///
/// Format:
/// ```text
/// c comment
/// p cnf <n_vars> <n_clauses>
/// 1 -2 3 0
/// -1 2 0
/// ```
pub fn parse_dimacs(input: &str) -> Result<Solver, String> {
    let mut n_vars: Option<u32> = None;
    let mut current: Vec<Lit> = Vec::new();
    let mut clauses: Vec<Vec<Lit>> = Vec::new();
    for (lineno, line) in input.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('c') {
            continue;
        }
        if line.starts_with("p ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 || parts[1] != "cnf" {
                return Err(format!("line {}: bad header `{line}`", lineno + 1));
            }
            n_vars = Some(
                parts[2]
                    .parse()
                    .map_err(|e| format!("line {}: {e}", lineno + 1))?,
            );
            continue;
        }
        // Literal stream; 0 terminates a clause.
        for tok in line.split_whitespace() {
            let l: Lit = tok
                .parse()
                .map_err(|e| format!("line {}: {e}", lineno + 1))?;
            if l == 0 {
                clauses.push(std::mem::take(&mut current));
            } else {
                current.push(l);
            }
        }
    }
    if !current.is_empty() {
        clauses.push(current);
    }
    let n = n_vars.ok_or("missing `p cnf` header")?;
    let mut solver = Solver::new(n);
    for c in clauses {
        if !solver.add_clause(c) {
            // Empty clause → UNSAT but still a valid instance.
            break;
        }
    }
    Ok(solver)
}

/// Parse an **extended DIMACS** string in which lines beginning with
/// `x` are XOR (parity) constraints, as emitted by CryptoMiniSat and by
/// the `EC-Index-Calculus-Benchmarks` generator.
///
/// The convention is that `x 1 2 3 0` means `x₁ ⊕ x₂ ⊕ x₃ = 1`, and
/// each negated literal flips the right-hand side, so
/// `x -1 2 3 0` means `x₁ ⊕ x₂ ⊕ x₃ = 0`.  Ordinary clause lines are
/// parsed exactly as in [`parse_dimacs`].
///
/// XOR rows go to [`Solver::add_xor`] rather than being Tseitin-
/// expanded, so the returned solver reasons about them by Gaussian
/// elimination.
pub fn parse_dimacs_xor(input: &str) -> Result<Solver, String> {
    let mut n_vars: Option<u32> = None;
    let mut clauses: Vec<Vec<Lit>> = Vec::new();
    let mut xors: Vec<(Vec<u32>, bool)> = Vec::new();

    for (lineno, line) in input.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('c') {
            continue;
        }
        if line.starts_with("p ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 || parts[1] != "cnf" {
                return Err(format!("line {}: bad header `{line}`", lineno + 1));
            }
            n_vars = Some(
                parts[2]
                    .parse()
                    .map_err(|e| format!("line {}: {e}", lineno + 1))?,
            );
            continue;
        }
        let is_xor = line.starts_with('x');
        let payload = if is_xor { &line[1..] } else { line };
        let mut lits: Vec<Lit> = Vec::new();
        for tok in payload.split_whitespace() {
            let l: Lit = tok
                .parse()
                .map_err(|_| format!("line {}: bad literal `{tok}`", lineno + 1))?;
            if l == 0 {
                break;
            }
            lits.push(l);
        }
        if lits.is_empty() {
            continue;
        }
        if is_xor {
            // rhs starts true and flips once per negated literal.
            let negations = lits.iter().filter(|l| **l < 0).count();
            let rhs = negations % 2 == 0;
            xors.push((lits.iter().map(|l| l.unsigned_abs()).collect(), rhs));
        } else {
            clauses.push(lits);
        }
    }

    let n = n_vars.ok_or("missing `p cnf` header")?;
    let mut solver = Solver::new(n);
    for c in clauses {
        if !solver.add_clause(c) {
            break;
        }
    }
    for (vars, rhs) in xors {
        if !solver.add_xor(&vars, rhs) {
            break;
        }
    }
    Ok(solver)
}

/// Emit a DIMACS string for a solver's *original* clauses.
pub fn to_dimacs(solver: &Solver) -> String {
    let mut s = format!("p cnf {} {}\n", solver.n_vars, solver.n_orig_clauses);
    for c in solver.clauses.iter().take(solver.n_orig_clauses) {
        for &l in c {
            s.push_str(&l.to_string());
            s.push(' ');
        }
        s.push_str("0\n");
    }
    s
}

/// Convenience: verify that a model satisfies all clauses.
pub fn check_model(clauses: &[Vec<Lit>], model: &[bool]) -> bool {
    for c in clauses {
        let mut sat = false;
        for &l in c {
            let v = var_of(l) as usize;
            let val = if is_neg(l) { !model[v] } else { model[v] };
            if val {
                sat = true;
                break;
            }
        }
        if !sat {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── native XOR reasoning ────────────────────────────────────────

    /// `x₁ ⊕ x₂ = 1` and `x₂ ⊕ x₃ = 1` and `x₁ ⊕ x₃ = 1` is the odd-
    /// cycle parity contradiction: summing all three rows gives `0 = 1`.
    /// Gaussian elimination must see it; resolution would have to search.
    #[test]
    fn xor_odd_cycle_is_unsat() {
        let mut s = Solver::new(3);
        assert!(s.add_xor(&[1, 2], true));
        assert!(s.add_xor(&[2, 3], true));
        assert!(s.add_xor(&[1, 3], true));
        assert_eq!(s.solve(), SolveResult::Unsat);
    }

    /// A triangular XOR system with a unique solution, derivable by
    /// propagation alone (no decisions needed).
    #[test]
    fn xor_system_propagates_to_unique_solution() {
        let mut s = Solver::new(3);
        s.add_xor(&[1], true); //           x₁ = 1
        s.add_xor(&[1, 2], false); //  x₁ ⊕ x₂ = 0  → x₂ = 1
        s.add_xor(&[2, 3], true); //   x₂ ⊕ x₃ = 1  → x₃ = 0
        assert_eq!(s.solve(), SolveResult::Sat);
        let m = s.model();
        assert_eq!((m[0], m[1], m[2]), (true, true, false));
        assert!(s.check_xors(&m));
    }

    /// Duplicated variables inside one row must cancel: `x ⊕ x ⊕ y = 1`
    /// is just `y = 1`.
    #[test]
    fn xor_duplicate_vars_cancel() {
        let mut s = Solver::new(2);
        assert!(s.add_xor(&[1, 1, 2], true));
        assert_eq!(s.solve(), SolveResult::Sat);
        assert!(s.model()[1], "y must be forced true");
    }

    /// `x ⊕ x = 1` reduces to `0 = 1` and is rejected at add time.
    #[test]
    fn xor_empty_row_with_true_rhs_is_unsat() {
        let mut s = Solver::new(2);
        assert!(!s.add_xor(&[1, 1], true));
        assert_eq!(s.solve(), SolveResult::Unsat);
    }

    /// XOR rows and CNF clauses must constrain each other: the parity
    /// rows admit two solutions, and the clause rules one out.
    #[test]
    fn xor_and_cnf_interact() {
        let mut s = Solver::new(3);
        s.add_xor(&[1, 2], false); // x₁ = x₂
        s.add_xor(&[2, 3], false); // x₂ = x₃
        s.add_clause(vec![1]); //     x₁ = 1  ⇒ all true
        assert_eq!(s.solve(), SolveResult::Sat);
        let m = s.model();
        assert_eq!((m[0], m[1], m[2]), (true, true, true));
        assert!(s.check_xors(&m));
    }

    /// **Differential test against brute force.**  Random dense parity
    /// systems, decided exhaustively and by the solver; every verdict
    /// must agree, and every SAT model must satisfy every row.
    #[test]
    fn xor_engine_agrees_with_brute_force() {
        let n_vars = 9u32;
        let mut state = 0x2545_F491_4F6C_DD1Du64; // xorshift64*
        let mut next = || {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545_F491_4F6C_DD1D)
        };

        for trial in 0..200 {
            let n_rows = 3 + (next() % 6) as usize;
            let mut rows: Vec<(Vec<u32>, bool)> = Vec::new();
            for _ in 0..n_rows {
                let mut vars: Vec<u32> = Vec::new();
                for v in 1..=n_vars {
                    if next() % 2 == 0 {
                        vars.push(v);
                    }
                }
                if vars.is_empty() {
                    vars.push(1 + (next() % n_vars as u64) as u32);
                }
                rows.push((vars, next() % 2 == 0));
            }

            // Ground truth by exhaustive search.
            let mut brute_sat = false;
            for a in 0..(1u32 << n_vars) {
                if rows.iter().all(|(vars, rhs)| {
                    let parity = vars.iter().filter(|v| (a >> (**v - 1)) & 1 == 1).count() % 2 == 1;
                    parity == *rhs
                }) {
                    brute_sat = true;
                    break;
                }
            }

            let mut s = Solver::new(n_vars);
            let mut trivially_unsat = false;
            for (vars, rhs) in &rows {
                if !s.add_xor(vars, *rhs) {
                    trivially_unsat = true;
                }
            }
            let res = s.solve();

            if brute_sat {
                assert_eq!(res, SolveResult::Sat, "trial {trial}: solver missed a model");
                let m = s.model();
                assert!(s.check_xors(&m), "trial {trial}: model violates a row");
            } else {
                assert_eq!(res, SolveResult::Unsat, "trial {trial}: solver invented a model");
                let _ = trivially_unsat;
            }
        }
    }

    /// The extended-DIMACS reader must apply CryptoMiniSat's parity
    /// convention: bare literals mean `= 1`, one negation flips to `= 0`.
    #[test]
    fn parse_dimacs_xor_reads_parity_convention() {
        // x₁ ⊕ x₂ ⊕ x₃ = 1 (no negations), x₁ ⊕ x₂ = 0 (one negation),
        // plus the unit clause x₃.  Together these force x₃ = 1 and
        // x₁ = x₂.
        let src = "p cnf 3 3\nx 1 2 3 0\nx -1 2 0\n3 0\n";
        let mut s = parse_dimacs_xor(src).expect("parse");
        assert_eq!(s.n_xors(), 2);
        assert_eq!(s.solve(), SolveResult::Sat);
        let m = s.model();
        assert!(s.check_xors(&m));
        assert!(m[2], "the unit clause forces x₃ true");
        assert!(m[0] ^ m[1] ^ m[2], "first row must have odd parity");
        assert!(!(m[0] ^ m[1]), "second row must have even parity");

        // Flipping the unit clause to ¬x₃ makes the same system UNSAT.
        let src_bad = "p cnf 3 3\nx 1 2 3 0\nx -1 2 0\n-3 0\n";
        assert_eq!(
            parse_dimacs_xor(src_bad).expect("parse").solve(),
            SolveResult::Unsat
        );
    }

    /// A small known-SAT instance.
    #[test]
    fn simple_sat() {
        let mut s = Solver::new(3);
        assert!(s.add_clause(vec![1, 2]));
        assert!(s.add_clause(vec![-1, 3]));
        assert!(s.add_clause(vec![-2, -3]));
        assert!(s.add_clause(vec![1, -3]));
        assert_eq!(s.solve(), SolveResult::Sat);
        let m = s.model();
        // Verify against the original clauses.
        let orig = vec![vec![1, 2], vec![-1, 3], vec![-2, -3], vec![1, -3]];
        assert!(check_model(&orig, &m));
    }

    /// Trivially UNSAT: x ∧ ¬x.  Note: `add_clause` may return false
    /// for the second unit clause because enqueuing `-1` conflicts
    /// with the already-propagated `1`; both that path and `solve()`
    /// must reach UNSAT.
    #[test]
    fn trivial_unsat() {
        let mut s = Solver::new(1);
        assert!(s.add_clause(vec![1]));
        // The contradictory unit may be detected immediately at add
        // time (returns false) or surface from solve(); either is
        // valid UNSAT detection.
        let added = s.add_clause(vec![-1]);
        if added {
            assert_eq!(s.solve(), SolveResult::Unsat);
        } else {
            // add_clause already detected UNSAT; solve() must agree.
            assert_eq!(s.solve(), SolveResult::Unsat);
        }
    }

    /// The pigeonhole principle PHP(3, 2): 3 pigeons in 2 holes is
    /// UNSAT. Variable `x_{i,j}` = pigeon `i` is in hole `j`. Indexed
    /// 1..=6 as: (i-1)*2 + j.
    #[test]
    fn pigeonhole_3_into_2_is_unsat() {
        let x = |i: u32, j: u32| ((i - 1) * 2 + j) as i32;
        let mut s = Solver::new(6);
        // Each pigeon in some hole.
        for i in 1..=3 {
            s.add_clause(vec![x(i, 1), x(i, 2)]);
        }
        // No two pigeons in the same hole.
        for j in 1..=2 {
            for i1 in 1..=3 {
                for i2 in (i1 + 1)..=3 {
                    s.add_clause(vec![-x(i1, j), -x(i2, j)]);
                }
            }
        }
        assert_eq!(s.solve(), SolveResult::Unsat);
    }

    /// Pigeonhole 5-into-4 is also UNSAT; non-trivial for CDCL.
    #[test]
    fn pigeonhole_5_into_4_is_unsat() {
        let n_pigeons = 5u32;
        let n_holes = 4u32;
        let x = |i: u32, j: u32| ((i - 1) * n_holes + j) as i32;
        let mut s = Solver::new(n_pigeons * n_holes);
        for i in 1..=n_pigeons {
            let mut c = Vec::with_capacity(n_holes as usize);
            for j in 1..=n_holes {
                c.push(x(i, j));
            }
            s.add_clause(c);
        }
        for j in 1..=n_holes {
            for i1 in 1..=n_pigeons {
                for i2 in (i1 + 1)..=n_pigeons {
                    s.add_clause(vec![-x(i1, j), -x(i2, j)]);
                }
            }
        }
        assert_eq!(s.solve(), SolveResult::Unsat);
    }

    /// Round-trip a CNF through DIMACS.
    #[test]
    fn dimacs_roundtrip() {
        let mut s = Solver::new(4);
        s.add_clause(vec![1, -2, 3]);
        s.add_clause(vec![-1, 2]);
        s.add_clause(vec![-3, 4]);
        let dimacs = to_dimacs(&s);
        assert!(dimacs.starts_with("p cnf 4 3"));
        let s2 = parse_dimacs(&dimacs).expect("parse failed");
        assert_eq!(s2.n_vars(), 4);
        assert_eq!(s2.n_clauses(), 3);
    }

    /// Luby sequence check.
    #[test]
    fn luby_sequence_starts_correctly() {
        let seq: Vec<u64> = (1..=15).map(luby).collect();
        assert_eq!(seq, vec![1, 1, 2, 1, 1, 2, 4, 1, 1, 2, 1, 1, 2, 4, 8]);
    }

    /// Random k-SAT instance at the satisfiability threshold (clause-
    /// to-variable ratio ≈ 4.26 for 3-SAT). Should mostly be SAT.
    #[test]
    fn random_3sat_around_threshold() {
        let n = 20u32;
        let m = (4.0 * n as f64) as usize;
        let mut state = 0xc0ffee_u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            state >> 33
        };
        let mut s = Solver::new(n);
        let mut orig_clauses = Vec::new();
        for _ in 0..m {
            let mut clause = Vec::with_capacity(3);
            let mut vars_in_clause = HashSet::new();
            while clause.len() < 3 {
                let v = ((next() as u32) % n) + 1;
                if vars_in_clause.insert(v) {
                    let sign = (next() & 1) as u32;
                    clause.push(if sign == 0 { v as i32 } else { -(v as i32) });
                }
            }
            orig_clauses.push(clause.clone());
            s.add_clause(clause);
        }
        let result = s.solve();
        if result == SolveResult::Sat {
            let m = s.model();
            assert!(check_model(&orig_clauses, &m));
        }
        // Either SAT or UNSAT is acceptable; what we're checking is
        // that the solver terminates.
        assert_ne!(result, SolveResult::Unknown);
    }
}
