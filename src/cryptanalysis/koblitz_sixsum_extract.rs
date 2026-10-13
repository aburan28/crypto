//! # Higher-arity (5- and 6-sum) relation extraction over the Frobenius pair table.
//!
//! The compact-orbit index calculus on a Koblitz curve finds **4-sum**
//! relations `R = ±σ^a(u) ± σ^b(u')` by scanning a table of pair sums
//! `u = P_i + P_j` (one entry per `⟨σ, −1⟩`-orbit of pairs, `K²n` entries
//! for `K` orbit columns): for every table state `s` (an entry and a
//! Frobenius shift) it solves the quadratic `S₃(x_R, x_s, X) = 0` for `X`
//! and looks the two roots up in the table by the canonical name of their
//! Frobenius orbit.  `RESEARCH_ECC2K130_IC_FEASIBILITY.md` §4.7 sketches a
//! **6-sum** extraction over the same table and asks what it would cost.
//! This module builds the three obvious kernels and counts what they do,
//! so that the question is settled by measurement rather than by a model:
//!
//! - **kernel A, "iterate one table"** ([`Kernel::IterateTable`]): for
//!   every signed, shifted state `s₁` form `R' = R − s₁` (one point
//!   addition) and run the 4-sum `S₃`-root probe of `R'`, stopping at the
//!   first verified hit — `R = s₁ + s₂ + s₃`, six base points;
//! - **kernel B, "S₄ per probe"** ([`Kernel::S4PerProbe`]): for every
//!   unordered pair of abscissa states `(s₁, s₂)` solve the quartic
//!   `S₄(x_R, x_{s₁}, x_{s₂}, X) = 0` over `F_{2^n}` (a real degree-4
//!   root-finding, [`s4_roots`]) and look its roots up in the table;
//! - **kernel C, the 5-sum control** ([`Kernel::FiveSumControl`]): `R − P_e`
//!   for every base point, then the 4-sum probe;
//! - and the **4-sum baseline** itself ([`Kernel::FourSum`]).
//!
//! Every candidate is lifted to points, its signs resolved by group
//! arithmetic, and the whole sum re-added in the single-word *and* in the
//! general curve arithmetic before it counts ([`verify_relation`]); a
//! candidate that does not close is a `lift_failures` tick, never a
//! relation.  The counters ([`ProbeCounters`]) keep the `S₃` solves, `S₄`
//! solves, point additions, table lookups and canonicalisations apart, so
//! a kernel can be compared in the ladder's own unit (`S₃` solves) and in
//! the feasibility note's unit (targets) at once.
//!
//! ## The unit question
//!
//! In `S₃`-solve units a 6-sum probe `(s₁, s₂)` of kernel A tests every
//! table state for `s₃` at once, exactly as a 4-sum probe `s₁` tests every
//! `s₂`: the hit probability per solve is the table's coverage of the
//! abscissa line, `≈ 4K²n²/r` per root, whatever the arity.  So in the
//! **dense regime** (several relations per target, first hit inside the
//! first target) the expected solves per relation are the same,
//! `r/(4A)` for `A = K²n²` abscissa states, for the 4-sum, the 5-sum
//! control and kernel A, and what the higher arity buys is only
//! *targets* per relation.  In the **sparse regime** (fewer than one
//! relation of that arity per target, `B^m/(m!·r) < 1`) a kernel pays
//! its complete scan on every empty target, so its cost is `μ ×` the
//! dense value with `μ` the number of times one relation is hit in a
//! complete scan (`μ = 6` for the 4-sum, 30 for the 5-sum, 90 for
//! kernel A, 45 for kernel B).  A 6-sum kernel is dense wherever
//! `B⁶/(720r) ≫ 1`, so against a sparse 4-sum it saves at most that
//! factor `6` in `S₃` solves, and nothing against a dense one.  The
//! note's `45·r/(n⁴K⁴)` counts first-summand probes assuming an `S₄`
//! oracle that tests all `(s₂, s₃)` pairs at once; no such oracle
//! exists, and kernel B's quartic — whose four roots are
//! `x(R ± s₁ ± s₂)` — is exactly the union of the two `S₃` root pairs
//! kernel A computes for `R ∓ s₁`, so the honest `S₄` solve is at best
//! half the solves at many times the cost each.
//! `research/ic_techniques_20261012/sixsum_extraction/` records the
//! measurement.
//!
//! Nothing here is a key-recovery, sub-rho or challenge-scale claim: the
//! fixtures are the public registered Koblitz rungs and the targets are
//! `[a]G` for a seeded `a`.

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

use serde::Serialize;

use crate::cryptanalysis::ic_boundary::{
    fast_map, BinaryGroup, BinaryInstance, CountedGroup, FactorBase, FastMap, GroupOps,
    OracleCounters,
};
use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};
use crate::cryptanalysis::koblitz_fast::{FastPoint, FrobeniusCanon};
use crate::cryptanalysis::koblitz_fast_arith::s3_x_roots;
use crate::cryptanalysis::semaev_decomp::Gf2;

// ── Kernels ────────────────────────────────────────────────────────

/// Which extraction loop runs over the table.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Kernel {
    /// The 4-sum baseline: one `S₃`-root probe of `R` over the table.
    FourSum,
    /// Kernel A: `R − s₁` for every signed shifted state, then the 4-sum
    /// probe.  Six base points.
    IterateTable,
    /// Kernel B: the quartic `S₄(x_R, x_{s₁}, x_{s₂}, X)` for every
    /// unordered pair of abscissa states, roots looked up.  Six points.
    /// Pairs are visited row by row (`s₁` fixed, `s₂` running), so one
    /// row is the 4-sum scan of `R ± s₁`.
    S4PerProbe,
    /// Kernel B with the pairs visited along a diagonal stride
    /// (`(i, j) → (i+1, j+1)`), so consecutive quartics share neither
    /// `s₁` nor `s₂`: the same solves as [`Kernel::S4PerProbe`] in a
    /// complete scan, but the hits of one target — clustered by row,
    /// since a row has hits only when `R ± s₁` has a 4-sum — are met in
    /// a de-clustered order.  Six points.  Possible only because the
    /// quartic needs abscissae alone (no `R − s₁` point addition).
    S4Strided,
    /// Kernel C: `R − P_e` for every base point, then the 4-sum probe.
    /// Five points.
    FiveSumControl,
}

impl Kernel {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "4" | "4sum" | "four" | "foursum" | "baseline" => Ok(Self::FourSum),
            "a" | "iterate" | "iterate_table" | "6a" => Ok(Self::IterateTable),
            "b" | "s4" | "s4_per_probe" | "6b" => Ok(Self::S4PerProbe),
            "b2" | "bs" | "s4_strided" | "strided" | "6b2" => Ok(Self::S4Strided),
            "c" | "5" | "5sum" | "five" | "control" => Ok(Self::FiveSumControl),
            other => Err(format!("unknown kernel `{other}` (4sum | a | b | b2 | c)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::FourSum => "4sum",
            Self::IterateTable => "6sum_a_iterate_table",
            Self::S4PerProbe => "6sum_b_s4_per_probe",
            Self::S4Strided => "6sum_b2_s4_strided",
            Self::FiveSumControl => "5sum_c_control",
        }
    }

    /// Base points in a relation from this kernel.
    pub fn arity(self) -> usize {
        match self {
            Self::FourSum => 4,
            Self::IterateTable | Self::S4PerProbe | Self::S4Strided => 6,
            Self::FiveSumControl => 5,
        }
    }
}

// ── Counters, budgets, relations ───────────────────────────────────

/// What a probe loop did, in the units the kernels are compared in.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct ProbeCounters {
    /// Targets handed to [`SixSumTable::extract`].
    pub targets: u64,
    /// Quadratic `S₃(x_R, x_s, X)` root solves ([`s3_x_roots`]).
    pub s3_solves: u64,
    /// Quartic `S₄(x_R, x₁, x₂, X)` root solves ([`s4_roots`]).
    pub s4_solves: u64,
    /// **Probes in the feasibility note's unit**: first-summand states
    /// tried.  For the 4-sum every abscissa state (one `S₃` solve each);
    /// for kernel A every signed shifted state `s₁` (each followed by a
    /// whole inner 4-sum scan); for kernel C every base point `P_e`; for
    /// kernel B every unordered pair `(s₁, s₂)` (one `S₄` solve each).
    /// §4.7's `45·r/(n⁴K⁴)` is a count in this unit.
    pub outer_probes: u64,
    /// Point additions on the probe side: `R − s₁`, state points, sign
    /// resolution.  Verification additions are counted apart.
    pub point_adds: u64,
    /// Canonical-orbit names computed for root lookups.
    pub canonicalisations: u64,
    /// Table lookups of a root's orbit name.
    pub lookups: u64,
    /// Lookups that found an entry (before sign resolution).
    pub candidate_hits: u64,
    /// Candidates whose signs could not be resolved into the target, or
    /// that failed the independent re-addition.  Fail-closed: these are
    /// never relations.
    pub lift_failures: u64,
    /// Relations that passed [`verify_relation`].
    pub verified_relations: u64,
    /// Relations rejected by [`verify_relation`] after sign resolution
    /// claimed they closed.  Must stay zero; a non-zero value is an
    /// arithmetic bug, not a yield.
    pub verification_rejects: u64,
    /// Additions spent re-adding verified relations (fast arithmetic).
    pub verification_adds: u64,
    /// `R − s₁ = O` or `R − P_e = O`: the target *is* a table state or a
    /// base point, a lower-arity relation this kernel does not emit.
    pub degenerate_skips: u64,
    /// Targets whose scan ran to completion (no early abort, no budget
    /// stop).
    pub complete_scans: u64,
    /// Targets stopped by the budget.
    pub budget_stops: u64,
    /// Wall inside [`SixSumTable::extract`].
    pub wall_ns: u64,
}

impl ProbeCounters {
    pub fn merge(&mut self, o: &ProbeCounters) {
        self.targets += o.targets;
        self.s3_solves += o.s3_solves;
        self.s4_solves += o.s4_solves;
        self.outer_probes += o.outer_probes;
        self.point_adds += o.point_adds;
        self.canonicalisations += o.canonicalisations;
        self.lookups += o.lookups;
        self.candidate_hits += o.candidate_hits;
        self.lift_failures += o.lift_failures;
        self.verified_relations += o.verified_relations;
        self.verification_rejects += o.verification_rejects;
        self.verification_adds += o.verification_adds;
        self.degenerate_skips += o.degenerate_skips;
        self.complete_scans += o.complete_scans;
        self.budget_stops += o.budget_stops;
        self.wall_ns += o.wall_ns;
    }
}

/// When a scan stops.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    /// Stop a target at its first verified relation.
    pub early_abort: bool,
    /// Cumulative `S₃` solves (over the counters handed in) after which
    /// the scan stops; `u64::MAX` for none.
    pub max_s3_solves: u64,
    /// The same for `S₄` solves.
    pub max_s4_solves: u64,
    /// Wall deadline for the scan.
    pub deadline: Option<Instant>,
}

impl Default for Budget {
    fn default() -> Self {
        Self {
            early_abort: true,
            max_s3_solves: u64::MAX,
            max_s4_solves: u64::MAX,
            deadline: None,
        }
    }
}

impl Budget {
    #[inline]
    fn exhausted(&self, ctr: &ProbeCounters) -> bool {
        ctr.s3_solves >= self.max_s3_solves
            || ctr.s4_solves >= self.max_s4_solves
            || self.deadline.is_some_and(|d| Instant::now() >= d)
    }
}

/// A verified relation `Σ P_{indices} = target`.
#[derive(Clone, Debug, Serialize)]
pub struct Relation {
    /// Base indices, sorted; the sum of these points is the target.
    pub indices: Vec<usize>,
    /// `indices.len()`.
    pub arity: usize,
    /// Distinct base points left after cancelling `P` against `−P`:
    /// the arity the relation really has as a group equation.
    pub support: usize,
    /// `(column, coefficient mod r)` of the relation row, zero
    /// coefficients dropped: `Σ coef·log(column) = log(target)` in the
    /// base's signed-Frobenius-orbit column model.
    pub row: Vec<(usize, u64)>,
    pub target_x: u64,
    pub target_y: u64,
}

/// Why a scan ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum ScanEnd {
    /// Every state was tried.
    Complete,
    /// The first verified relation stopped it (early abort).
    FirstHit,
    /// The budget stopped it.
    Budget,
}

// ── The table ──────────────────────────────────────────────────────

/// One pair-sum entry, `u = P_a + P_b`, with the canonical name of its
/// abscissa's Frobenius orbit computed at the shift `shift`.
#[derive(Clone, Copy, Debug)]
pub struct PairState {
    pub point: FastPoint,
    pub a: u32,
    pub b: u32,
    /// `t_u` with `canon = rotl^{t_u}(coords(x_u))`.
    pub shift: u32,
}

/// What the build did.
#[derive(Clone, Debug, Default, Serialize)]
pub struct BuildStats {
    pub base_points: u64,
    pub columns: u64,
    pub orbit_representatives: u64,
    /// Pairs `(rep, j)` added.
    pub pair_adds: u64,
    /// Entries kept: one per distinct orbit name of a pair-sum abscissa.
    pub states: u64,
    /// Pairs whose sum's orbit name an earlier pair already had.
    pub duplicate_keys: u64,
    pub canonicalisations: u64,
    pub frobenius_maps: u64,
    /// `[coef_i]·rep(col_i) == P_i` checks of the column model.
    pub coefficient_checks: u64,
    pub wall_ns: u64,
}

/// **The Frobenius pair table with its states materialised**: one entry
/// per `⟨σ, −1⟩`-orbit of pairs, the sum point kept so a state can be
/// used as a summand (`R − s₁`) and not only as a key.  Same orbit and
/// representative logic as `ic_boundary::FrobeniusPairTable`; the
/// abscissae of every entry's `n` Frobenius shifts are tabulated so a
/// probe reads `x(σ^t u)` instead of squaring.
pub struct SixSumTable {
    pub states: Vec<PairState>,
    /// `shifted_x[e·n + t] = x(σ^t u_e)`.
    shifted_x: Vec<u64>,
    /// Orbit name of `x_u` → state index.
    map: FastMap<u32>,
    canon: FrobeniusCanon,
    /// `member_of[i·n + t]`: the base index of `σ^t(P_i)`.
    member_of: Vec<u32>,
    n: u32,
    b: u64,
    pub build: BuildStats,
}

impl SixSumTable {
    /// Build over a base closed under the signed Frobenius (the
    /// `compact-orbit-scan` base).  Checks the base's column model —
    /// every point is its column representative times its coefficient —
    /// so that a verified point relation is also a verified row.
    pub fn build(inst: &BinaryInstance, fb: &FactorBase<FastPoint>) -> Result<Self, String> {
        let start = Instant::now();
        let n = inst.n as usize;
        let f = fb.points.len();
        if n == 0 || f == 0 {
            return Err("empty base or field".into());
        }
        if f >= u32::MAX as usize / n.max(1) {
            return Err("base too wide for u32 indices".into());
        }
        let canon = FrobeniusCanon::new(&inst.gf, inst.n)
            .ok_or("no normal basis found for the Frobenius canonicalisation")?;
        let mut stats = BuildStats {
            base_points: f as u64,
            columns: fb.columns as u64,
            ..BuildStats::default()
        };
        // σ^t(P_i) by index.
        let mut member_of = vec![0u32; f * n];
        for i in 0..f {
            let mut p = fb.points[i];
            member_of[i * n] = i as u32;
            for t in 1..n {
                p = inst.fast.frobenius_k(p, 1);
                stats.frobenius_maps += 1;
                let idx = *fb
                    .point_index
                    .get(&p.pack())
                    .ok_or_else(|| format!("base is not closed under the Frobenius at point {i}"))?;
                member_of[i * n + t] = idx as u32;
            }
        }
        // The column model: P_i = [coef_i / coef_rep] · P_rep for the first
        // point of its column.
        let group = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let mut rep_of_col: Vec<Option<usize>> = vec![None; fb.columns];
        for i in 0..f {
            let c = fb.col_of[i];
            if rep_of_col[c].is_none() {
                rep_of_col[c] = Some(i);
            }
        }
        for i in 0..f {
            let rep = rep_of_col[fb.col_of[i]].expect("every column has its first point");
            let lhs = group.mul(&mut ops, fb.points[rep], fb.coef_of[i]);
            let rhs = group.mul(&mut ops, fb.points[i], fb.coef_of[rep]);
            stats.coefficient_checks += 1;
            if lhs != rhs {
                return Err(format!(
                    "column model fails at point {i}: [coef]·rep != [coef_rep]·P"
                ));
            }
        }
        // Orbits of base points by canonical abscissa name; one
        // representative each.
        let mut orbit_of = vec![0u32; f];
        let mut rep_of_orbit: Vec<usize> = Vec::new();
        {
            let mut orbit_index: FastMap<u32> = fast_map(f / n + 1);
            for i in 0..f {
                let key = canon.canon(fb.points[i].x);
                stats.canonicalisations += 1;
                let next = rep_of_orbit.len() as u32;
                let id = *orbit_index.entry(key).or_insert(next);
                if id == next {
                    rep_of_orbit.push(i);
                }
                orbit_of[i] = id;
            }
        }
        stats.orbit_representatives = rep_of_orbit.len() as u64;
        let mut states: Vec<PairState> = Vec::with_capacity(f * f / (4 * n) + 1);
        let mut map: FastMap<u32> = fast_map(f * f / (4 * n) + 1);
        for (orbit, &rep) in rep_of_orbit.iter().enumerate() {
            for j in 0..f {
                if (orbit_of[j] as usize) < orbit {
                    continue;
                }
                let u = inst.fast.add(fb.points[rep], fb.points[j]);
                stats.pair_adds += 1;
                if u.infinity {
                    continue;
                }
                let (key, t_u) = canon.canon_with_shift(u.x);
                stats.canonicalisations += 1;
                let next = states.len() as u32;
                let id = *map.entry(key).or_insert(next);
                if id == next {
                    states.push(PairState {
                        point: u,
                        a: rep as u32,
                        b: j as u32,
                        shift: t_u,
                    });
                } else {
                    stats.duplicate_keys += 1;
                }
            }
        }
        stats.states = states.len() as u64;
        let mut shifted_x = vec![0u64; states.len() * n];
        for (e, s) in states.iter().enumerate() {
            let mut x = s.point.x;
            for t in 0..n {
                shifted_x[e * n + t] = x;
                x = inst.gf.sqr(x);
            }
        }
        stats.wall_ns = start.elapsed().as_nanos() as u64;
        Ok(Self {
            states,
            shifted_x,
            map,
            canon,
            member_of,
            n: inst.n,
            b: inst.b,
            build: stats,
        })
    }

    /// Entries in the table (`≈ K²n`).
    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// Abscissa states, `states × n` shifts.
    pub fn abscissa_states(&self) -> usize {
        self.states.len() * self.n as usize
    }

    /// Signed shifted states, `states × n × 2`: the summands of kernel A.
    pub fn signed_states(&self) -> usize {
        self.abscissa_states() * 2
    }

    #[inline]
    fn x_of(&self, e: usize, t: usize) -> u64 {
        self.shifted_x[e * self.n as usize + t]
    }

    /// The state `(e, t)` whose shifted abscissa is `x`, if the table has
    /// its orbit: `x = x(σ^t u_e)`.
    #[inline]
    fn lookup(&self, x: u64, ctr: &mut ProbeCounters) -> Option<(usize, usize)> {
        ctr.canonicalisations += 1;
        let (key, t_x) = self.canon.canon_with_shift(x);
        ctr.lookups += 1;
        let &e = self.map.get(&key)?;
        ctr.candidate_hits += 1;
        let n = self.n;
        let t_u = self.states[e as usize].shift;
        // canon = rotl^{t_u}(coords(x_u)) = rotl^{t_x}(coords(x)), the
        // Frobenius being rotl by one, so x = σ^t(x_u), t = t_u − t_x.
        let t = ((t_u + n - t_x) % n) as usize;
        Some((e as usize, t))
    }

    /// The base indices of `σ^t(u_e) = P_a' + P_b'`.
    #[inline]
    fn state_indices(&self, e: usize, t: usize) -> [usize; 2] {
        let n = self.n as usize;
        let s = &self.states[e];
        [
            self.member_of[s.a as usize * n + t] as usize,
            self.member_of[s.b as usize * n + t] as usize,
        ]
    }

    /// `σ^t(u_e)` as a point, by adding its two members (one addition;
    /// exact, no squaring chain).
    #[inline]
    fn state_point(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        e: usize,
        t: usize,
        ctr: &mut ProbeCounters,
    ) -> (FastPoint, [usize; 2]) {
        let idx = self.state_indices(e, t);
        ctr.point_adds += 1;
        (inst.fast.add(fb.points[idx[0]], fb.points[idx[1]]), idx)
    }

    /// Given `target = ±s ± s'` for the state points `s = σ^{t}(u_e)` and
    /// `s' = σ^{t'}(u_{e'})`, the four signed base indices, or `None`
    /// when no sign choice closes.
    fn resolve_pair_signs(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        target: FastPoint,
        (e1, t1): (usize, usize),
        (e2, t2): (usize, usize),
        ctr: &mut ProbeCounters,
    ) -> Option<[usize; 4]> {
        let (s1, i1) = self.state_point(inst, fb, e1, t1, ctr);
        let (s2, i2) = self.state_point(inst, fb, e2, t2, ctr);
        let neg = |idx: [usize; 2]| [fb.neg_index[idx[0]], fb.neg_index[idx[1]]];
        ctr.point_adds += 1;
        let q = inst.fast.add(target, inst.fast.neg(s1));
        if q == s2 {
            return Some([i1[0], i1[1], i2[0], i2[1]]);
        }
        if q == inst.fast.neg(s2) {
            let m = neg(i2);
            return Some([i1[0], i1[1], m[0], m[1]]);
        }
        ctr.point_adds += 1;
        let q = inst.fast.add(target, s1);
        let m1 = neg(i1);
        if q == s2 {
            return Some([m1[0], m1[1], i2[0], i2[1]]);
        }
        if q == inst.fast.neg(s2) {
            let m2 = neg(i2);
            return Some([m1[0], m1[1], m2[0], m2[1]]);
        }
        None
    }

    /// **The 4-sum `S₃`-root probe** of `target` over the table, from
    /// state `start` cyclically.  For every abscissa state `(e, t)` it
    /// solves `S₃(x_R, x(σ^t u_e), X) = 0`, looks both roots up, resolves
    /// the signs of a hit and hands the four base indices to `sink`,
    /// which returns `false` to stop.  Returns how the scan ended.  The
    /// sink's indices are *not* yet verified by re-addition; callers go
    /// through [`verify_relation`].
    pub fn probe_four_sum(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        target: FastPoint,
        start: usize,
        budget: &Budget,
        ctr: &mut ProbeCounters,
        sink: &mut dyn FnMut(&[usize], &mut ProbeCounters) -> bool,
    ) -> ScanEnd {
        if target.infinity || self.states.is_empty() {
            return ScanEnd::Complete;
        }
        let gf = &inst.gf;
        let n = self.n as usize;
        let count = self.states.len();
        let start = start % count;
        for step in 0..count {
            let e = (start + step) % count;
            if budget.exhausted(ctr) {
                return ScanEnd::Budget;
            }
            for t in 0..n {
                let xs = self.x_of(e, t);
                ctr.s3_solves += 1;
                let Some(roots) = s3_x_roots(gf, self.b, target.x, xs) else {
                    continue;
                };
                for (k, &root) in roots.iter().enumerate() {
                    if k == 1 && roots[1] == roots[0] {
                        break;
                    }
                    let Some((e2, t2)) = self.lookup(root, ctr) else {
                        continue;
                    };
                    match self.resolve_pair_signs(inst, fb, target, (e, t), (e2, t2), ctr) {
                        Some(idx) => {
                            if !sink(&idx, ctr) {
                                return ScanEnd::FirstHit;
                            }
                        }
                        None => ctr.lift_failures += 1,
                    }
                }
            }
        }
        ScanEnd::Complete
    }

    /// **Run one kernel on one target.**  `start` rotates the scan order
    /// (the example derives it from the seed and the target index, so
    /// relations do not all share the first state); every relation
    /// returned has passed [`verify_relation`].  With `budget.early_abort`
    /// the scan stops at the first verified relation.
    pub fn extract(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        target: FastPoint,
        kernel: Kernel,
        start: u64,
        budget: &Budget,
        ctr: &mut ProbeCounters,
    ) -> (Vec<Relation>, ScanEnd) {
        let begin = Instant::now();
        ctr.targets += 1;
        let mut out: Vec<Relation> = Vec::new();
        let mut seen: HashMap<Vec<usize>, ()> = HashMap::new();
        let early = budget.early_abort;
        // Verify, dedupe, record; `false` stops the scan.
        let mut accept = |indices: &[usize], ctr: &mut ProbeCounters, out: &mut Vec<Relation>| -> bool {
            let mut sorted = indices.to_vec();
            sorted.sort_unstable();
            if seen.contains_key(&sorted) {
                return true;
            }
            match verify_relation(inst, fb, &sorted, target, ctr) {
                Ok(rel) => {
                    seen.insert(sorted, ());
                    ctr.verified_relations += 1;
                    out.push(rel);
                    !early
                }
                Err(_) => {
                    ctr.verification_rejects += 1;
                    true
                }
            }
        };
        let end = match kernel {
            Kernel::FourSum => {
                let count = self.states.len().max(1);
                let s3_before = ctr.s3_solves;
                let end = self.probe_four_sum(
                    inst,
                    fb,
                    target,
                    (start % count as u64) as usize,
                    budget,
                    ctr,
                    &mut |idx, c| accept(idx, c, &mut out),
                );
                ctr.outer_probes += ctr.s3_solves - s3_before;
                end
            }
            Kernel::IterateTable => {
                let n = self.n as usize;
                let total = self.signed_states();
                if total == 0 || target.infinity {
                    ScanEnd::Complete
                } else {
                    let start = (start % total as u64) as usize;
                    let mut end = ScanEnd::Complete;
                    for step in 0..total {
                        if budget.exhausted(ctr) {
                            end = ScanEnd::Budget;
                            break;
                        }
                        let s = (start + step) % total;
                        let sign = s & 1;
                        let et = s >> 1;
                        let (e1, t1) = (et / n, et % n);
                        let (mut s1, mut i1) = self.state_point(inst, fb, e1, t1, ctr);
                        if sign == 1 {
                            s1 = inst.fast.neg(s1);
                            i1 = [fb.neg_index[i1[0]], fb.neg_index[i1[1]]];
                        }
                        ctr.outer_probes += 1;
                        ctr.point_adds += 1;
                        let rest = inst.fast.add(target, inst.fast.neg(s1));
                        if rest.infinity {
                            ctr.degenerate_skips += 1;
                            continue;
                        }
                        let inner_start = (e1 + 1) % self.states.len().max(1);
                        let inner = self.probe_four_sum(
                            inst,
                            fb,
                            rest,
                            inner_start,
                            budget,
                            ctr,
                            &mut |idx, c| {
                                let full = [i1[0], i1[1], idx[0], idx[1], idx[2], idx[3]];
                                accept(&full, c, &mut out)
                            },
                        );
                        match inner {
                            ScanEnd::FirstHit => {
                                end = ScanEnd::FirstHit;
                                break;
                            }
                            ScanEnd::Budget => {
                                end = ScanEnd::Budget;
                                break;
                            }
                            ScanEnd::Complete => {}
                        }
                    }
                    end
                }
            }
            Kernel::FiveSumControl => {
                let total = fb.points.len();
                if total == 0 || target.infinity {
                    ScanEnd::Complete
                } else {
                    let start = (start % total as u64) as usize;
                    let mut end = ScanEnd::Complete;
                    for step in 0..total {
                        if budget.exhausted(ctr) {
                            end = ScanEnd::Budget;
                            break;
                        }
                        let i = (start + step) % total;
                        ctr.outer_probes += 1;
                        ctr.point_adds += 1;
                        let rest = inst.fast.add(target, inst.fast.neg(fb.points[i]));
                        if rest.infinity {
                            ctr.degenerate_skips += 1;
                            continue;
                        }
                        let inner_start = step % self.states.len().max(1);
                        let inner = self.probe_four_sum(
                            inst,
                            fb,
                            rest,
                            inner_start,
                            budget,
                            ctr,
                            &mut |idx, c| {
                                let full = [i, idx[0], idx[1], idx[2], idx[3]];
                                accept(&full, c, &mut out)
                            },
                        );
                        match inner {
                            ScanEnd::FirstHit => {
                                end = ScanEnd::FirstHit;
                                break;
                            }
                            ScanEnd::Budget => {
                                end = ScanEnd::Budget;
                                break;
                            }
                            ScanEnd::Complete => {}
                        }
                    }
                    end
                }
            }
            Kernel::S4PerProbe => {
                self.probe_s4_pairs(inst, fb, target, start, false, budget, ctr, &mut |idx, c| {
                    accept(idx, c, &mut out)
                })
            }
            Kernel::S4Strided => {
                self.probe_s4_pairs(inst, fb, target, start, true, budget, ctr, &mut |idx, c| {
                    accept(idx, c, &mut out)
                })
            }
        };
        match end {
            ScanEnd::Complete => ctr.complete_scans += 1,
            ScanEnd::Budget => ctr.budget_stops += 1,
            ScanEnd::FirstHit => {}
        }
        ctr.wall_ns += begin.elapsed().as_nanos() as u64;
        (out, end)
    }

    /// Kernel B's loop: every unordered pair of abscissa states
    /// `(i ≤ j)` over the flattened `(e, t)` index; one quartic solve
    /// per pair; each root looked up and the eight sign choices of
    /// `R = ±s₁ ± s₂ ± s₃` resolved by group arithmetic.  Row order
    /// (`strided = false`): outer `i` rotated by `start`, `j` running
    /// `i..total`.  Diagonal order (`strided = true`): the cells of the
    /// `total × total` square visited as `idx = start + k·(total + 1)
    /// mod total²`, i.e. `(i, j) → (i + 1, j + 1)`, cells with `j < i`
    /// skipped without a solve; `total + 1` is coprime to `total²`, so a
    /// complete walk visits every pair once.
    #[allow(clippy::too_many_arguments)]
    fn probe_s4_pairs(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        target: FastPoint,
        start: u64,
        strided: bool,
        budget: &Budget,
        ctr: &mut ProbeCounters,
        sink: &mut dyn FnMut(&[usize], &mut ProbeCounters) -> bool,
    ) -> ScanEnd {
        let total = self.abscissa_states();
        if total == 0 || target.infinity {
            return ScanEnd::Complete;
        }
        let gf = &inst.gf;
        let n = self.n as usize;
        let mut roots: Vec<u64> = Vec::with_capacity(4);
        // One quartic for the pair `(i, j)`; `Some(false)` stops.
        let solve_pair = |i: usize,
                              j: usize,
                              roots: &mut Vec<u64>,
                              ctr: &mut ProbeCounters,
                              sink: &mut dyn FnMut(&[usize], &mut ProbeCounters) -> bool|
         -> bool {
            let (e1, t1) = (i / n, i % n);
            let (e2, t2) = (j / n, j % n);
            let x1 = self.x_of(e1, t1);
            let x2 = self.x_of(e2, t2);
            ctr.s4_solves += 1;
            ctr.outer_probes += 1;
            let quartic = s4_in_x4(gf, self.b, target.x, x1, x2);
            roots.clear();
            s4_roots_split(gf, &quartic, roots);
            roots.sort_unstable();
            roots.dedup();
            for &root in roots.iter() {
                let Some((e3, t3)) = self.lookup(root, ctr) else {
                    continue;
                };
                match self.resolve_triple_signs(inst, fb, target, (e1, t1), (e2, t2), (e3, t3), ctr)
                {
                    Some(idx) => {
                        if !sink(&idx, ctr) {
                            return false;
                        }
                    }
                    None => ctr.lift_failures += 1,
                }
            }
            true
        };
        if strided {
            let square = (total as u128) * (total as u128);
            let stride = total as u128 + 1;
            let mut idx = (start as u128) % square;
            let mut k: u128 = 0;
            while k < square {
                let i = (idx / total as u128) as usize;
                let j = (idx % total as u128) as usize;
                idx += stride;
                if idx >= square {
                    idx -= square;
                }
                k += 1;
                if j < i {
                    continue;
                }
                if budget.exhausted(ctr) {
                    return ScanEnd::Budget;
                }
                if !solve_pair(i, j, &mut roots, ctr, sink) {
                    return ScanEnd::FirstHit;
                }
            }
        } else {
            let start = (start % total as u64) as usize;
            for step in 0..total {
                let i = (start + step) % total;
                for j in i..total {
                    if budget.exhausted(ctr) {
                        return ScanEnd::Budget;
                    }
                    if !solve_pair(i, j, &mut roots, ctr, sink) {
                        return ScanEnd::FirstHit;
                    }
                }
            }
        }
        ScanEnd::Complete
    }

    /// `target = ±s₁ ± s₂ ± s₃` → six signed base indices, or `None`.
    #[allow(clippy::too_many_arguments)]
    fn resolve_triple_signs(
        &self,
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        target: FastPoint,
        a: (usize, usize),
        b: (usize, usize),
        c: (usize, usize),
        ctr: &mut ProbeCounters,
    ) -> Option<[usize; 6]> {
        let (s1, i1) = self.state_point(inst, fb, a.0, a.1, ctr);
        let (s2, i2) = self.state_point(inst, fb, b.0, b.1, ctr);
        let (s3, i3) = self.state_point(inst, fb, c.0, c.1, ctr);
        let neg = |idx: [usize; 2]| [fb.neg_index[idx[0]], fb.neg_index[idx[1]]];
        for sign1 in 0..2 {
            let (p1, j1) = if sign1 == 0 {
                (s1, i1)
            } else {
                (inst.fast.neg(s1), neg(i1))
            };
            ctr.point_adds += 1;
            let q1 = inst.fast.add(target, inst.fast.neg(p1));
            for sign2 in 0..2 {
                let (p2, j2) = if sign2 == 0 {
                    (s2, i2)
                } else {
                    (inst.fast.neg(s2), neg(i2))
                };
                ctr.point_adds += 1;
                let q = inst.fast.add(q1, inst.fast.neg(p2));
                if q == s3 {
                    return Some([j1[0], j1[1], j2[0], j2[1], i3[0], i3[1]]);
                }
                if q == inst.fast.neg(s3) {
                    let m = neg(i3);
                    return Some([j1[0], j1[1], j2[0], j2[1], m[0], m[1]]);
                }
            }
        }
        None
    }
}

// ── Verification ───────────────────────────────────────────────────

/// **Fail-closed check of a candidate relation.**  The base points named
/// by `indices` are re-added in the single-word arithmetic and, when the
/// instance carries the general `KoblitzCurve`, independently in the
/// general arithmetic; both sums must equal `target`, every point must
/// be on the curve, and the target must be a finite point.  Returns the
/// relation with its row over the base's columns and its support after
/// cancellation.
pub fn verify_relation(
    inst: &BinaryInstance,
    fb: &FactorBase<FastPoint>,
    indices: &[usize],
    target: FastPoint,
    ctr: &mut ProbeCounters,
) -> Result<Relation, String> {
    if target.infinity {
        return Err("target is the identity".into());
    }
    if indices.is_empty() {
        return Err("empty relation".into());
    }
    let mut sum = FastPoint::INFINITY;
    for &i in indices {
        let p = *fb
            .points
            .get(i)
            .ok_or_else(|| format!("index {i} outside the base"))?;
        if !inst.fast.is_on_curve(p) {
            return Err(format!("base point {i} is off the curve"));
        }
        ctr.verification_adds += 1;
        sum = inst.fast.add(sum, p);
    }
    if sum != target {
        return Err("fast re-addition does not reach the target".into());
    }
    if let Some(kc) = inst.koblitz.as_ref() {
        let mut general = crate::binary_ecc::BinaryPoint::Infinity;
        for &i in indices {
            general = kc.add(&general, &inst.fast.lower(fb.points[i]));
        }
        if general != inst.fast.lower(target) {
            return Err("general-arithmetic re-addition does not reach the target".into());
        }
    }
    // Support after cancelling P against −P.
    let mut counts: HashMap<u64, i64> = HashMap::new();
    for &i in indices {
        let p = fb.points[i];
        let key = p.pack();
        let neg_key = inst.fast.neg(p).pack();
        match counts.get_mut(&neg_key) {
            Some(c) if *c > 0 => *c -= 1,
            _ => *counts.entry(key).or_insert(0) += 1,
        }
    }
    let support = counts.values().map(|&c| c.max(0) as usize).sum();
    let mut row: BTreeMap<usize, u64> = BTreeMap::new();
    for &i in indices {
        let c = fb.col_of[i];
        let v = row.entry(c).or_insert(0);
        *v = ((*v as u128 + fb.coef_of[i] as u128) % inst.r as u128) as u64;
    }
    let row: Vec<(usize, u64)> = row.into_iter().filter(|&(_, v)| v != 0).collect();
    Ok(Relation {
        indices: indices.to_vec(),
        arity: indices.len(),
        support,
        row,
        target_x: target.x,
        target_y: target.y,
    })
}

// ── The quartic S₄ and its roots ───────────────────────────────────

/// Coefficients `c[0..=4]` of `S₄(x₁, x₂, x₃, X)` as a polynomial in `X`
/// over `F_{2^n}` for `y² + xy = x³ + ax² + b`:
///
/// ```text
/// S₄ = Res_Y(S₃(x₁, x₂, Y), S₃(x₃, X, Y))
///    = (A₁C₂ + A₂C₁)² + (A₁B₂ + A₂B₁)(B₁C₂ + B₂C₁)
/// ```
///
/// with `A₁ = (x₁+x₂)², B₁ = x₁x₂, C₁ = (x₁x₂)² + b` constants and
/// `A₂ = x₃² + X², B₂ = x₃X, C₂ = x₃²X² + b` polynomials in `X`
/// (`binary_semaev::binary_semaev_s4`'s formula, expanded in its last
/// argument).  Degree 4 in general, lower on degenerate inputs.
pub fn s4_in_x4(gf: &Gf2, b: u64, x1: u64, x2: u64, x3: u64) -> [u64; 5] {
    let b = b & gf.mask;
    let a1 = gf.sqr(x1 ^ x2);
    let b1 = gf.mul(x1, x2);
    let c1 = gf.sqr(b1) ^ b;
    let x3s = gf.sqr(x3);
    // T = A₁C₂ + A₂C₁ = [A₁b + C₁x₃², 0, A₁x₃² + C₁]
    let t0 = gf.mul(a1, b) ^ gf.mul(c1, x3s);
    let t2 = gf.mul(a1, x3s) ^ c1;
    // U = A₁B₂ + A₂B₁ = [B₁x₃², A₁x₃, B₁]
    let u0 = gf.mul(b1, x3s);
    let u1 = gf.mul(a1, x3);
    let u2 = b1;
    // V = B₁C₂ + B₂C₁ = [B₁b, C₁x₃, B₁x₃²]
    let v0 = gf.mul(b1, b);
    let v1 = gf.mul(c1, x3);
    let v2 = gf.mul(b1, x3s);
    // T² (char 2: coefficients squared, degrees doubled) + U·V.
    let mut c = [0u64; 5];
    c[0] = gf.sqr(t0);
    c[4] = gf.sqr(t2);
    c[0] ^= gf.mul(u0, v0);
    c[1] ^= gf.mul(u0, v1) ^ gf.mul(u1, v0);
    c[2] ^= gf.mul(u0, v2) ^ gf.mul(u1, v1) ^ gf.mul(u2, v0);
    c[3] ^= gf.mul(u1, v2) ^ gf.mul(u2, v1);
    c[4] ^= gf.mul(u2, v2);
    c
}

/// Evaluate a polynomial (low coefficient first) at `x`.
pub fn poly_eval(gf: &Gf2, coeffs: &[u64], x: u64) -> u64 {
    let mut acc = 0u64;
    for &c in coeffs.iter().rev() {
        acc = gf.mul(acc, x) ^ c;
    }
    acc
}

fn trim(p: &mut Vec<u64>) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

fn degree(p: &[u64]) -> Option<usize> {
    p.len().checked_sub(1)
}

/// Make monic; `p` non-zero.
fn monic(gf: &Gf2, p: &mut [u64]) {
    let lead = *p.last().expect("non-zero polynomial");
    if lead != 1 {
        let inv = gf.inv(lead);
        for c in p.iter_mut() {
            *c = gf.mul(*c, inv);
        }
    }
}

/// `a mod f`, `f` monic and non-constant.
fn rem(gf: &Gf2, a: &[u64], f: &[u64]) -> Vec<u64> {
    let df = degree(f).expect("non-zero modulus");
    let mut r = a.to_vec();
    trim(&mut r);
    while let Some(dr) = degree(&r) {
        if dr < df {
            break;
        }
        let coef = r[dr];
        let shift = dr - df;
        for (k, &fk) in f.iter().enumerate() {
            r[shift + k] ^= gf.mul(coef, fk);
        }
        trim(&mut r);
    }
    r
}

/// `(a div f, a mod f)`, `f` monic and non-constant.
fn divmod(gf: &Gf2, a: &[u64], f: &[u64]) -> (Vec<u64>, Vec<u64>) {
    let df = degree(f).expect("non-zero modulus");
    let mut r = a.to_vec();
    trim(&mut r);
    let mut q = vec![0u64; r.len().saturating_sub(df) + 1];
    while let Some(dr) = degree(&r) {
        if dr < df {
            break;
        }
        let coef = r[dr];
        let shift = dr - df;
        q[shift] ^= coef;
        for (k, &fk) in f.iter().enumerate() {
            r[shift + k] ^= gf.mul(coef, fk);
        }
        trim(&mut r);
    }
    trim(&mut q);
    (q, r)
}

/// `p² mod f` in characteristic 2.
fn square_mod(gf: &Gf2, p: &[u64], f: &[u64]) -> Vec<u64> {
    if p.is_empty() {
        return Vec::new();
    }
    let mut sq = vec![0u64; 2 * p.len() - 1];
    for (i, &c) in p.iter().enumerate() {
        sq[2 * i] = gf.sqr(c);
    }
    rem(gf, &sq, f)
}

/// Monic `gcd(a, b)`; the zero polynomial is empty.
fn gcd(gf: &Gf2, a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    trim(&mut a);
    trim(&mut b);
    while !b.is_empty() {
        monic(gf, &mut b);
        let r = if b.len() == 1 { Vec::new() } else { rem(gf, &a, &b) };
        a = b;
        b = r;
    }
    if !a.is_empty() {
        monic(gf, &mut a);
    }
    a
}

/// Split a monic polynomial known to be a product of distinct linear
/// factors over `F_{2^n}` into its roots, by the trace: `Tr(βX) mod g`
/// for the basis elements `β = 2^i` in turn until a proper factor
/// appears, then recurse (deterministic Berlekamp trace splitting).
fn split_linear(gf: &Gf2, g: &[u64], roots: &mut Vec<u64>) {
    let Some(d) = degree(g) else {
        return;
    };
    match d {
        0 => {}
        1 => roots.push(g[0]),
        _ => {
            let n = gf.n;
            for i in 0..n {
                let beta = 1u64 << i;
                // acc = Σ_{j<n} (βX)^{2^j} mod g
                let mut term = rem(gf, &[0, beta], g);
                let mut acc = term.clone();
                for _ in 1..n {
                    term = square_mod(gf, &term, g);
                    if acc.len() < term.len() {
                        acc.resize(term.len(), 0);
                    }
                    for (a, &t) in acc.iter_mut().zip(&term) {
                        *a ^= t;
                    }
                    trim(&mut acc);
                }
                let h = gcd(gf, g, &acc);
                if let Some(dh) = degree(&h) {
                    if dh > 0 && dh < d {
                        let (q, _) = divmod(gf, g, &h);
                        split_linear(gf, &h, roots);
                        split_linear(gf, &q, roots);
                        return;
                    }
                }
            }
            // No basis element separates the roots: `g` is not a product
            // of distinct linear factors over the field.  Report nothing.
        }
    }
}

/// **Every distinct root in `F_{2^n}`** of the polynomial `coeffs` (low
/// coefficient first): `gcd(f, X^{2^n} − X)` by `n` squarings modulo `f`,
/// then trace splitting.  The zero polynomial and constants have no
/// roots reported.
pub fn s4_roots(gf: &Gf2, coeffs: &[u64], roots: &mut Vec<u64>) {
    let mut f = coeffs.to_vec();
    trim(&mut f);
    let Some(d) = degree(&f) else {
        return;
    };
    if d == 0 {
        return;
    }
    monic(gf, &mut f);
    if d == 1 {
        roots.push(f[0]);
        return;
    }
    // h = X^{2^n} mod f
    let mut h = rem(gf, &[0, 1], &f);
    for _ in 0..gf.n {
        h = square_mod(gf, &h, &f);
    }
    // h + X
    if h.len() < 2 {
        h.resize(2, 0);
    }
    h[1] ^= 1;
    trim(&mut h);
    let g = gcd(gf, &f, &h);
    split_linear(gf, &g, roots);
}

/// The roots of `coeffs` **assuming it splits completely over the
/// field** (as `S₄(x_R, x₁, x₂, X)` does when `x_R, x₁, x₂` are abscissae
/// of rational points: its roots are `x(R ± s₁ ± s₂)`): trace splitting
/// without the `X^{2^n} − X` gcd, quadratics finished by the half-trace.
/// Every root reported is checked by evaluation, so on a polynomial that
/// does not split the result is a (possibly empty) subset of the true
/// roots, never a wrong one.  This is what kernel B runs; [`s4_roots`]
/// is the general routine the tests compare it with.
///
/// Cost on a quartic: `n` squarings modulo `f` for the residues
/// `X^{2^j} mod f` (reused by every trace), then per basis element `β`
/// one trace `Σ β^{2^j} X^{2^j}` and one small gcd, then two half-trace
/// quadratics — a dozen or so `S₃`-solve equivalents.
pub fn s4_roots_split(gf: &Gf2, coeffs: &[u64], roots: &mut Vec<u64>) {
    let mut f = coeffs.to_vec();
    trim(&mut f);
    let Some(d) = degree(&f) else {
        return;
    };
    if d == 0 {
        return;
    }
    monic(gf, &mut f);
    let before = roots.len();
    if d == 1 {
        roots.push(f[0]);
    } else if d == 2 {
        quadratic_roots(gf, &f, roots);
    } else {
        // Squarefree part first: repeated roots defeat the trace split.
        // In characteristic 2 the derivative keeps the odd-degree terms.
        let mut deriv = vec![0u64; f.len() - 1];
        for (i, &c) in f.iter().enumerate().skip(1) {
            if i % 2 == 1 {
                deriv[i - 1] = c;
            }
        }
        trim(&mut deriv);
        let g = if deriv.is_empty() {
            // f = h²: the roots of h.
            let mut h = vec![0u64; f.len().div_ceil(2)];
            for (i, slot) in h.iter_mut().enumerate() {
                *slot = gf.sqr_k(f[2 * i], gf.n - 1);
            }
            trim(&mut h);
            h
        } else {
            let common = gcd(gf, &f, &deriv);
            if degree(&common).unwrap_or(0) > 0 {
                divmod(gf, &f, &common).0
            } else {
                f.clone()
            }
        };
        match degree(&g) {
            Some(1) => roots.push(g[0]),
            Some(2) => quadratic_roots(gf, &g, roots),
            Some(dg) if dg >= 3 => {
                // Residues X^{2^j} mod g, j = 0..n.
                let n = gf.n as usize;
                let mut powers: Vec<Vec<u64>> = Vec::with_capacity(n);
                let mut p = rem(gf, &[0, 1], &g);
                for _ in 0..n {
                    powers.push(p.clone());
                    p = square_mod(gf, &p, &g);
                }
                split_with_powers(gf, &g, &powers, 0, roots);
            }
            _ => {}
        }
    }
    // Keep only genuine roots of the polynomial handed in.
    let mut k = before;
    while k < roots.len() {
        if poly_eval(gf, coeffs, roots[k]) != 0 {
            roots.swap_remove(k);
        } else {
            k += 1;
        }
    }
}

/// Roots in `F_{2^n}` of a monic quadratic `X² + aX + b`: `X = a·z` with
/// `z² + z = b/a²` (Artin–Schreier, one half-trace for odd `n`); `a = 0`
/// gives the double root `√b`.  Roots are checked by evaluation.
fn quadratic_roots(gf: &Gf2, g: &[u64], roots: &mut Vec<u64>) {
    debug_assert_eq!(g.len(), 3);
    let (b, a) = (g[0], g[1]);
    if a == 0 {
        roots.push(gf.sqr_k(b, gf.n - 1));
        return;
    }
    let inv_a = gf.inv(a);
    let c = gf.mul(b, gf.sqr(inv_a));
    let Some(z) = crate::cryptanalysis::koblitz_fast_arith::artin_schreier_root(gf, gf.n, c)
    else {
        return;
    };
    let r1 = gf.mul(a, z);
    let r2 = r1 ^ a;
    for r in [r1, r2] {
        if poly_eval(gf, g, r) == 0 {
            roots.push(r);
        }
    }
}

/// Trace splitting of a monic squarefree `g` of degree `≥ 3` whose
/// residues `powers[j] = X^{2^j} mod g` are in hand: for the basis
/// elements `β = 2^i`, `i ≥ from`, form `Tr(βX) mod g = Σ β^{2^j}·powers[j]`
/// and split by `gcd(g, Tr(βX))`; the factors inherit the residues by
/// reduction and recurse, degree-2 factors going to the half-trace.
fn split_with_powers(gf: &Gf2, g: &[u64], powers: &[Vec<u64>], from: u32, roots: &mut Vec<u64>) {
    let Some(d) = degree(g) else {
        return;
    };
    match d {
        0 => return,
        1 => {
            roots.push(g[0]);
            return;
        }
        2 => {
            quadratic_roots(gf, g, roots);
            return;
        }
        _ => {}
    }
    let n = gf.n;
    for i in from..n {
        let mut beta = 1u64 << i;
        let mut trace = vec![0u64; d];
        for pw in powers {
            for (slot, &c) in trace.iter_mut().zip(pw) {
                *slot ^= gf.mul(beta, c);
            }
            beta = gf.sqr(beta);
        }
        trim(&mut trace);
        let h = gcd(gf, g, &trace);
        if let Some(dh) = degree(&h) {
            if dh > 0 && dh < d {
                let (q, _) = divmod(gf, g, &h);
                let reduce = |m: &[u64]| -> Vec<Vec<u64>> {
                    powers.iter().map(|p| rem(gf, p, m)).collect()
                };
                let ph = reduce(&h);
                split_with_powers(gf, &h, &ph, i + 1, roots);
                let pq = reduce(&q);
                split_with_powers(gf, &q, &pq, i + 1, roots);
                return;
            }
        }
    }
}

// ── Framework oracle ───────────────────────────────────────────────

/// The kernels as a framework [`DecompositionOracle`]: `prepare` builds
/// the table, `decompose` returns the first verified relation's indices
/// (early abort).  Parameters: `kernel=4sum|a|b|c`.
pub struct SixSumOracle<'i> {
    pub kernel: Kernel,
    pub instance: &'i BinaryInstance,
    table: Option<SixSumTable>,
    pub counters: ProbeCounters,
    budget: Budget,
    calls: u64,
}

impl<'i> SixSumOracle<'i> {
    pub fn new(kernel: Kernel, instance: &'i BinaryInstance) -> Self {
        Self {
            kernel,
            instance,
            table: None,
            counters: ProbeCounters::default(),
            budget: Budget::default(),
            calls: 0,
        }
    }

    /// Cap the `S₃`/`S₄` solves one `decompose` call may spend.
    pub fn with_per_call_budget(mut self, max_s3: u64, max_s4: u64) -> Self {
        self.budget.max_s3_solves = max_s3;
        self.budget.max_s4_solves = max_s4;
        self
    }

    pub fn table(&self) -> Option<&SixSumTable> {
        self.table.as_ref()
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for SixSumOracle<'_> {
    fn name(&self) -> &str {
        "sixsum-extract"
    }

    fn summands(&self) -> u32 {
        self.kernel.arity() as u32
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "{} over the materialised Frobenius pair table: S3-root probes with every candidate re-added before it counts",
            self.kernel.name()
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[(
            "kernel",
            "4sum | a (iterate table) | b (S4 per probe, row order) | b2 (S4 per probe, diagonal stride) | c (5-sum control)",
        )]
    }

    fn prepare(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        if let Some(k) = params.get("kernel") {
            self.kernel = Kernel::parse(k)?;
        }
        let table = SixSumTable::build(self.instance, fb)?;
        ops.adds += table.build.pair_adds;
        self.table = Some(table);
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut counts = BTreeMap::new();
        if let Some(t) = &self.table {
            counts.insert("pair_table_states".into(), t.build.states);
            counts.insert("canonicalisations".into(), t.build.canonicalisations);
            counts.insert("frobenius_maps".into(), t.build.frobenius_maps);
            counts.insert("coefficient_checks".into(), t.build.coefficient_checks);
        }
        counts
    }

    fn decompose(
        &mut self,
        _ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        self.calls += 1;
        let mut budget = self.budget;
        budget.early_abort = true;
        // Per-call budgets are relative to this call's own counters.
        let mut local = ProbeCounters::default();
        let (rels, _) = table.extract(
            self.instance,
            fb,
            point,
            self.kernel,
            self.calls,
            &budget,
            &mut local,
        );
        ops.adds += local.point_adds + local.verification_adds;
        counters.lookups += local.lookups;
        counters.canonicalisations += local.canonicalisations;
        counters.lift_failures += local.lift_failures;
        counters.sqrt_solves += local.s3_solves;
        counters.s4_pairs += local.s4_solves;
        self.counters.merge(&local);
        rels.into_iter().next().map(|r| r.indices)
    }
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;
    use crate::cryptanalysis::ic_framework::plugins::CompactOrbitScanBase;
    use crate::cryptanalysis::ic_framework::stages::FactorBaseBuilder;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn base(inst: &BinaryInstance, columns: usize) -> FactorBase<FastPoint> {
        let group = BinaryGroup(&inst.fast);
        let ctx = InstanceCtx {
            group: &group,
            generator: inst.generator,
            target: group.identity(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(inst.n),
        };
        let mut params = Params::default();
        params.set("columns", columns.to_string());
        params.set("raw_x_cap", "100000000");
        let mut ops = GroupOps::default();
        CompactOrbitScanBase { instance: inst }
            .build(&ctx, &params, &mut ops)
            .expect("compact-orbit base")
    }

    fn sum_of(inst: &BinaryInstance, fb: &FactorBase<FastPoint>, idx: &[usize]) -> FastPoint {
        idx.iter()
            .fold(FastPoint::INFINITY, |s, &i| inst.fast.add(s, fb.points[i]))
    }

    fn brute_roots(gf: &Gf2, coeffs: &[u64]) -> Vec<u64> {
        let mut out: Vec<u64> = (0..(1u64 << gf.n))
            .filter(|&x| poly_eval(gf, coeffs, x) == 0)
            .collect();
        out.sort_unstable();
        out
    }

    /// The quartic agrees with the crate's reference `S₄` evaluated in
    /// the general arithmetic, and both root finders agree with brute
    /// force over the whole field at `n = 17`.
    #[test]
    fn s4_quartic_matches_reference_and_roots_match_brute_force() {
        use crate::cryptanalysis::binary_semaev::binary_semaev_s4;
        let inst = koblitz_instance(1, 17).expect("n=17 Koblitz instance");
        let gf = &inst.gf;
        let irr = &inst.irreducible;
        let b_el = gf.to_element(inst.b);
        let mut rng = StdRng::seed_from_u64(17);
        for round in 0..12 {
            let (x1, x2, x3) = (
                rng.gen::<u64>() & gf.mask,
                rng.gen::<u64>() & gf.mask,
                rng.gen::<u64>() & gf.mask,
            );
            let q = s4_in_x4(gf, inst.b, x1, x2, x3);
            for _ in 0..8 {
                let x = rng.gen::<u64>() & gf.mask;
                let reference = binary_semaev_s4(
                    &gf.to_element(x1),
                    &gf.to_element(x2),
                    &gf.to_element(x3),
                    &gf.to_element(x),
                    &b_el,
                    irr,
                );
                assert_eq!(
                    gf.from_element(&reference),
                    poly_eval(gf, &q, x),
                    "round {round}: S4 evaluation disagrees with the reference"
                );
            }
            let brute = brute_roots(gf, &q);
            let mut general = Vec::new();
            s4_roots(gf, &q, &mut general);
            general.sort_unstable();
            assert_eq!(general, brute, "round {round}: gcd root finder vs brute force");
            let mut split = Vec::new();
            s4_roots_split(gf, &q, &mut split);
            split.sort_unstable();
            // The split path reports only genuine roots and, when the
            // quartic splits, all of them.
            assert!(split.iter().all(|r| brute.contains(r)));
            if brute.len() == 4 {
                assert_eq!(split, brute, "round {round}: split path on a split quartic");
            }
        }
        // Polynomials built from rational points split completely and
        // the roots are exactly x(R ± s1 ± s2).
        let g = inst.generator;
        for k in 2..14u64 {
            let r = inst.fast.mul_u64(g, k);
            let s1 = inst.fast.mul_u64(g, 3 * k + 1);
            let s2 = inst.fast.mul_u64(g, 7 * k + 5);
            let q = s4_in_x4(gf, inst.b, r.x, s1.x, s2.x);
            let mut want: Vec<u64> = [
                inst.fast.add(inst.fast.add(r, s1), s2),
                inst.fast.add(inst.fast.add(r, s1), inst.fast.neg(s2)),
                inst.fast.add(inst.fast.add(r, inst.fast.neg(s1)), s2),
                inst.fast.add(inst.fast.add(r, inst.fast.neg(s1)), inst.fast.neg(s2)),
            ]
            .iter()
            .filter(|p| !p.infinity)
            .map(|p| p.x)
            .collect();
            want.sort_unstable();
            want.dedup();
            let mut split = Vec::new();
            s4_roots_split(gf, &q, &mut split);
            split.sort_unstable();
            assert_eq!(split, want, "k = {k}: roots are x(R ± s1 ± s2)");
            let mut general = Vec::new();
            s4_roots(gf, &q, &mut general);
            general.sort_unstable();
            assert_eq!(general, want, "k = {k}: general root finder");
        }
    }

    /// The table's states are its distinct keys, every state point is
    /// the sum of its two members, every member-shift is a base point,
    /// and the lookup of a shifted state abscissa returns that state and
    /// shift.
    #[test]
    fn the_table_states_and_lookups_are_consistent() {
        let inst = koblitz_instance(0, 23).expect("n=23 K_0 instance");
        let fb = base(&inst, 3);
        let table = SixSumTable::build(&inst, &fb).expect("table");
        assert_eq!(table.len() as u64, table.build.states);
        assert_eq!(table.map.len(), table.len());
        assert_eq!(fb.points.len(), 2 * 23 * 3);
        let mut ctr = ProbeCounters::default();
        for (e, s) in table.states.iter().enumerate() {
            assert_eq!(
                inst.fast.add(fb.points[s.a as usize], fb.points[s.b as usize]),
                s.point
            );
            for t in 0..23usize {
                let (p, _) = table.state_point(&inst, &fb, e, t, &mut ctr);
                assert_eq!(p, inst.fast.frobenius_k(s.point, t as u32), "e={e} t={t}");
                assert_eq!(p.x, table.x_of(e, t));
                let (e2, t2) = table.lookup(p.x, &mut ctr).expect("state abscissa is in the table");
                // Another entry may share the orbit (duplicate keys are
                // dropped), but the shifted abscissa must agree.
                assert_eq!(table.x_of(e2, t2), p.x, "e={e} t={t}");
            }
        }
        assert_eq!(ctr.lift_failures, 0);
    }

    /// Planted 4-, 5- and 6-sum targets built from the table's own
    /// states are recovered by the 4-sum probe, kernel C, kernel A and
    /// kernel B; the recovered relation is a true group equation.
    #[test]
    fn planted_targets_are_recovered_by_every_kernel() {
        let inst = koblitz_instance(0, 41).expect("n=41 K_0 instance");
        let fb = base(&inst, 2);
        let table = SixSumTable::build(&inst, &fb).expect("table");
        let n = 41usize;
        let mut ctr = ProbeCounters::default();
        // States used as summands: (e, t, sign).
        let pt = |e: usize, t: usize, neg: bool| -> (FastPoint, [usize; 2]) {
            let idx = table.state_indices(e, t);
            let p = inst.fast.add(fb.points[idx[0]], fb.points[idx[1]]);
            if neg {
                (inst.fast.neg(p), [fb.neg_index[idx[0]], fb.neg_index[idx[1]]])
            } else {
                (p, idx)
            }
        };
        let (s1, i1) = pt(0, 0, false);
        let (s2, i2) = pt(5 % table.len(), 7, true);
        let (s3, i3) = pt(11 % table.len(), 30, false);
        let budget = Budget::default();

        // 4-sum: R = s2 + s3.
        let r4 = inst.fast.add(s2, s3);
        let (rels, end) = table.extract(&inst, &fb, r4, Kernel::FourSum, 0, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit);
        assert_eq!(rels.len(), 1);
        assert_eq!(rels[0].arity, 4);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r4);

        // Kernel A: R = s1 + s2 + s3 with s1 the first state in scan order.
        let r6 = inst.fast.add(s1, r4);
        let (rels, end) = table.extract(&inst, &fb, r6, Kernel::IterateTable, 0, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit, "kernel A");
        assert_eq!(rels.len(), 1);
        assert_eq!(rels[0].arity, 6);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r6);
        let mut planted: Vec<usize> = [i1, i2, i3].concat();
        planted.sort_unstable();
        assert_eq!(rels[0].indices, planted, "kernel A recovers the planted six points");

        // Kernel B on the same target: pair (s1, s2) is reached from the
        // outer index 0.
        let (rels, end) = table.extract(&inst, &fb, r6, Kernel::S4PerProbe, 0, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit, "kernel B");
        assert_eq!(rels.len(), 1);
        assert_eq!(rels[0].arity, 6);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r6);
        assert_eq!(rels[0].indices, planted, "kernel B recovers the planted six points");
        assert!(ctr.s4_solves > 0);

        // Kernel B2 (diagonal stride) on the same target, started at
        // the cell (s1, s2) itself: the first quartic is the planted
        // pair.
        let a_states = table.abscissa_states() as u64;
        let j2 = (5 % table.len() as u64) * 41 + 7;
        let (rels, end) = table.extract(&inst, &fb, r6, Kernel::S4Strided, j2, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit, "kernel B2");
        assert_eq!(rels.len(), 1);
        assert_eq!(rels[0].arity, 6);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r6);
        assert_eq!(rels[0].indices, planted, "kernel B2 recovers the planted six points");
        // Started five diagonal steps before the cell (s2, s3), the walk
        // reaches it on its sixth quartic.
        let j3 = (11 % table.len() as u64) * 41 + 30;
        assert!(j3 > j2 + 5);
        let before = ctr.s4_solves;
        let start_cell = (j2 - 5) * a_states + (j3 - 5);
        let (rels, end) =
            table.extract(&inst, &fb, r6, Kernel::S4Strided, start_cell, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit, "kernel B2 along the diagonal");
        assert_eq!(ctr.s4_solves - before, 6, "six quartics along the diagonal");
        assert_eq!(rels.len(), 1);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r6);
        assert_eq!(rels[0].indices, planted);

        // Kernel C: R = P_0 + s2 + s3.
        let r5 = inst.fast.add(fb.points[0], r4);
        let (rels, end) = table.extract(&inst, &fb, r5, Kernel::FiveSumControl, 0, &budget, &mut ctr);
        assert_eq!(end, ScanEnd::FirstHit, "kernel C");
        assert_eq!(rels.len(), 1);
        assert_eq!(rels[0].arity, 5);
        assert_eq!(sum_of(&inst, &fb, &rels[0].indices), r5);
        let _ = n;
        assert_eq!(ctr.verification_rejects, 0);
        assert!(ctr.verified_relations >= 5);
        // Every row coefficient is a column sum of the base's coefficients.
        for rel in [&rels[0]] {
            for &(c, v) in &rel.row {
                assert!(c < fb.columns);
                assert!(v < inst.r);
            }
        }
    }

    /// A target whose logarithm is outside the 4-, 5- and 6-sum scalar
    /// lattice of a one-column base (checked independently on the
    /// coefficients) yields no relation from a complete scan of the
    /// 4-sum, 5-sum and kernel-A loops and none from a budgeted kernel-B
    /// scan; a target off the prime-order subgroup yields none either.
    #[test]
    fn non_decomposable_targets_yield_no_relation() {
        let inst = koblitz_instance(0, 41).expect("n=41 K_0 instance");
        let fb = base(&inst, 1);
        let table = SixSumTable::build(&inst, &fb).expect("table");
        let r = inst.r;
        let add = |a: u64, b: u64| ((a as u128 + b as u128) % r as u128) as u64;
        let sub = |a: u64, b: u64| ((a as u128 + r as u128 - b as u128) % r as u128) as u64;
        // C1 = coefficients, C2 = C1 + C1.
        let c1: std::collections::HashSet<u64> = fb.coef_of.iter().copied().collect();
        let mut c2: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for &a in &fb.coef_of {
            for &b in &fb.coef_of {
                c2.insert(add(a, b));
            }
        }
        let c2v: Vec<u64> = c2.iter().copied().collect();
        let in_lattice = |c: u64| -> bool {
            // 4-sum: c ∈ C2 + C2; 5-sum: C2 + C2 + C1; 6-sum: C2 + C2 + C2.
            for &a in &c2v {
                let rest = sub(c, a);
                if c2.contains(&rest) {
                    return true;
                }
                for &b in &c2v {
                    let rest2 = sub(rest, b);
                    if c2.contains(&rest2) || c1.contains(&rest2) {
                        return true;
                    }
                }
            }
            false
        };
        let rep = fb.points[0];
        let mut rng = StdRng::seed_from_u64(41);
        let mut c = rng.gen_range(2..r);
        let mut tries = 0;
        while in_lattice(c) {
            c = rng.gen_range(2..r);
            tries += 1;
            assert!(tries < 50, "could not find a scalar outside the lattice");
        }
        let target = inst.fast.mul_u64(rep, c);
        let mut ctr = ProbeCounters::default();
        let complete = Budget {
            early_abort: false,
            ..Budget::default()
        };
        for kernel in [Kernel::FourSum, Kernel::FiveSumControl, Kernel::IterateTable] {
            let (rels, end) = table.extract(&inst, &fb, target, kernel, 3, &complete, &mut ctr);
            assert!(rels.is_empty(), "{kernel:?} found a relation outside the lattice");
            assert_eq!(end, ScanEnd::Complete, "{kernel:?}");
        }
        let budgeted = Budget {
            early_abort: false,
            max_s4_solves: ctr.s4_solves + 20_000,
            ..Budget::default()
        };
        let (rels, _) = table.extract(&inst, &fb, target, Kernel::S4PerProbe, 3, &budgeted, &mut ctr);
        assert!(rels.is_empty(), "kernel B found a relation outside the lattice");
        let budgeted = Budget {
            early_abort: false,
            max_s4_solves: ctr.s4_solves + 20_000,
            ..Budget::default()
        };
        let (rels, _) = table.extract(&inst, &fb, target, Kernel::S4Strided, 3, &budgeted, &mut ctr);
        assert!(rels.is_empty(), "kernel B2 found a relation outside the lattice");
        assert_eq!(ctr.verified_relations, 0);
        assert_eq!(ctr.verification_rejects, 0);

        // Off the prime-order subgroup: R + T for the 2-torsion point T =
        // (0, √b).  No sum of subgroup points reaches it; every candidate
        // must fail sign resolution or verification.
        let torsion = FastPoint::affine(0, inst.gf.sqr_k(inst.b, inst.n - 1));
        assert!(inst.fast.is_on_curve(torsion));
        let off = inst.fast.add(inst.fast.mul_u64(rep, 12345), torsion);
        assert!(!inst.fast.mul_u64(off, r).infinity, "R + T is outside the subgroup");
        let mut ctr2 = ProbeCounters::default();
        for kernel in [Kernel::FourSum, Kernel::FiveSumControl, Kernel::IterateTable] {
            let (rels, _) = table.extract(&inst, &fb, off, kernel, 0, &complete, &mut ctr2);
            assert!(rels.is_empty(), "{kernel:?} verified a relation off the subgroup");
        }
        assert_eq!(ctr2.verified_relations, 0);
    }

    /// Verification is fail-closed: a wrong target or a tampered index
    /// list is rejected.
    #[test]
    fn verification_rejects_what_does_not_close() {
        let inst = koblitz_instance(0, 23).expect("n=23 K_0 instance");
        let fb = base(&inst, 2);
        let idx = [0usize, 3, 5, 8];
        let target = sum_of(&inst, &fb, &idx);
        let mut ctr = ProbeCounters::default();
        let rel = verify_relation(&inst, &fb, &idx, target, &mut ctr).expect("true relation");
        assert_eq!(rel.arity, 4);
        assert!(verify_relation(&inst, &fb, &idx, inst.fast.add(target, fb.points[1]), &mut ctr).is_err());
        assert!(verify_relation(&inst, &fb, &[0, 3, 5, 9], target, &mut ctr).is_err());
        assert!(verify_relation(&inst, &fb, &[0, 3, 5, usize::MAX], target, &mut ctr).is_err());
        // P + (−P) cancels in the support.
        let cancelling = [0usize, fb.neg_index[0], 3, 5];
        let t2 = sum_of(&inst, &fb, &cancelling);
        let rel = verify_relation(&inst, &fb, &cancelling, t2, &mut ctr).expect("relation");
        assert_eq!(rel.arity, 4);
        assert_eq!(rel.support, 2);
    }

    /// Random targets at a saturated small rung: the 4-sum baseline and
    /// kernel A both find verified relations and no verification is ever
    /// rejected.
    #[test]
    fn random_targets_yield_verified_relations_at_a_small_rung() {
        let inst = koblitz_instance(0, 23).expect("n=23 K_0 instance");
        let fb = base(&inst, 5);
        let table = SixSumTable::build(&inst, &fb).expect("table");
        let mut rng = StdRng::seed_from_u64(5);
        let mut ctr = ProbeCounters::default();
        let budget = Budget::default();
        let mut found = [0u32; 2];
        for k in 0..6u64 {
            let target = inst.fast.mul_u64(inst.generator, rng.gen_range(1..inst.r));
            for (slot, kernel) in [Kernel::FourSum, Kernel::IterateTable].into_iter().enumerate() {
                let (rels, _) = table.extract(&inst, &fb, target, kernel, k, &budget, &mut ctr);
                for rel in &rels {
                    assert_eq!(sum_of(&inst, &fb, &rel.indices), target);
                    assert_eq!(rel.arity, kernel.arity());
                    found[slot] += 1;
                }
            }
        }
        assert!(found[0] > 0, "4-sum found nothing at a saturated rung");
        assert!(found[1] > 0, "kernel A found nothing at a saturated rung");
        assert_eq!(ctr.verification_rejects, 0);
    }
}
