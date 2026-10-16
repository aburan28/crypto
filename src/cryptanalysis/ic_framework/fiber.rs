//! # Fiber-aware relation generation: the cofactor fiber of a target
//!
//! Every pipeline here writes its relations over the cofactor projection
//! `[h]P` of raw factor-base points (`collect_and_solve_with_completion`
//! multiplies the target coefficients by `h`), so a subgroup target `T`
//! stands for its whole **cofactor fiber** `{T + K : K ∈ E(F_q)[h]}`:
//! a decomposition `Σ P_i = T + K` is as good a row as `Σ P_i = T`,
//! because `[h]K = O`.  The oracles that landed before this module are
//! *fiber-blind*: they accept the one lift `T` and discard the other
//! `h − 1`.  This module is the fiber-aware family of
//! `research/isogeny_conductor_gap_ic_20261010/PROTOCOL.md` part A:
//!
//! | mode | what the oracle does per target |
//! |:--|:--|
//! | `lifts` | probes (or descends) every lift `T + K`: `h` oracle calls |
//! | `closed` | probes once against a pair table that holds `P_i + P_j + K` for every `K`: `h` times the entries, one probe |
//! | `combined` | (algebraic) one Weil-descended system with the target abscissa free and the fiber polynomial adjoined |
//!
//! The rational `h`-torsion `H = E(F_q)[h]` is computed in `prepare`
//! from `[r]P` over base points and charged to oracle set-up; when
//! `gcd(h, r) = 1` it has exactly `h` elements, and the oracle reports
//! how many it found (`fiber_size`) so a short fiber is visible rather
//! than silent.  With `h = 1` every mode is the blind oracle.
//!
//! Every hit is checked by the oracle itself — the recovered points are
//! added and compared with the lift they claim — so the relation loop
//! sees exact decompositions as it does from the blind oracles.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Duration;

use super::plugins::{charge_lift, charge_system_build, lift_probes, SolverHarness};
use super::stages::{
    BooleanSystem, DecompositionOracle, InstanceCtx, Params, SolverTotals, SolverVerdict,
    SystemShape, SystemSolver,
};
use crate::cryptanalysis::ic_boundary::{
    lift_abscissae, BinaryGroup, BinaryInstance, CountedGroup, FactorBase, FrobeniusPairTable,
    GroupOps, OracleCounters, PairTable,
};
use crate::cryptanalysis::koblitz_fast::{FastPoint, FrobeniusCanon};
use crate::cryptanalysis::pq_descent_symbolic::{
    descend, descend_fiber, max_n_prime, SymbolicDescent,
};

// ── the fiber ──────────────────────────────────────────────────────

/// `E(F_q)[h]`, the identity first, with how it was found.
#[derive(Clone, Debug)]
pub struct CofactorFiber<E> {
    pub points: Vec<E>,
    /// The cofactor the instance declares.
    pub h: u64,
    /// `points.len() == h`: the whole fiber was generated.
    pub complete: bool,
    /// Base points whose `[r]P` was taken.
    pub samples: u64,
    /// Probes and inserts of the generation's key set, for `bitop/v1`
    /// (`ecbench::work`) only.
    pub key_probes: u64,
}

/// The fiber generation's key-set probes, as table insert attempts.
fn add_fiber_key_probes<E>(c: &mut BTreeMap<String, u64>, fiber: Option<&CofactorFiber<E>>) {
    if let Some(f) = fiber {
        *c.entry("table_insert_attempts".to_string()).or_insert(0) += f.key_probes;
    }
}

impl<E> CofactorFiber<E> {
    pub fn size(&self) -> usize {
        self.points.len()
    }
}

/// Generate `E(F_q)[h]` from `[r]P` over the base points, charging the
/// scalar multiplications and closure additions to `ops`.  Stops at `h`
/// elements or when every base point has been used.
pub fn cofactor_fiber<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    r: u64,
    h: u64,
    ops: &mut GroupOps,
) -> CofactorFiber<G::Elt> {
    let mut points = vec![g.identity()];
    let mut keys: HashSet<u64> = HashSet::new();
    keys.insert(g.key(&g.identity()));
    let mut key_probes = 1u64;
    let mut samples = 0u64;
    if h > 1 {
        for &p in &fb.points {
            if points.len() as u64 >= h {
                break;
            }
            samples += 1;
            let k = g.mul(ops, p, r);
            key_probes += 1;
            if keys.contains(&g.key(&k)) {
                continue;
            }
            // Close the group under the new element: the cosets
            // `old + jK`, `j = 1 .. ord(K) − 1`.
            let old: Vec<G::Elt> = points.clone();
            let mut mult = k;
            loop {
                if g.is_identity(&mult) {
                    break;
                }
                for &o in &old {
                    let q = g.add(ops, o, mult);
                    key_probes += 1;
                    if keys.insert(g.key(&q)) {
                        points.push(q);
                    }
                }
                mult = g.add(ops, mult, k);
                if points.len() as u64 > h {
                    break;
                }
            }
        }
    }
    CofactorFiber {
        complete: points.len() as u64 == h,
        points,
        h,
        samples,
        key_probes,
    }
}

/// How the fiber is generated: `fiber_gen=legacy|coset`, `fiber_scale=d`.
///
/// `legacy` (the default) is [`cofactor_fiber`] exactly, so every committed
/// session replays with its counts.  Its closure walks every multiple
/// `jK`, `j < ord(K)`, of a new generator and adds each to every old
/// element, even once `jK` has fallen back into the old subgroup, where
/// no addition can produce a new point: `ord(K)·|old|` additions where
/// `[⟨old, K⟩ : old]·|old|` suffice.  On a non-cyclic `E(F_q)[h]` that
/// takes several generators, the excess is of the order of `h²/2`.
/// `coset` stops at the first multiple that lies in the old subgroup.
///
/// `fiber_scale=d` takes the subgroup `[d]E(F_q)[h]`, generated by
/// `[r·d]P` over base points, instead of the whole fiber: a fiber-closed
/// table of `|[d]E(F_q)[h]|` entries per pair, for the partial-fiber
/// sweep of round 2.  Its size is not known in advance, so every base point
/// is sampled (`complete` then means "every base point was used").
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FiberBuild {
    pub coset: bool,
    pub scale: u64,
}

impl Default for FiberBuild {
    fn default() -> Self {
        Self {
            coset: false,
            scale: 1,
        }
    }
}

impl FiberBuild {
    pub fn parse(params: &Params) -> Result<Self, String> {
        let coset = match params.get("fiber_gen").unwrap_or("legacy") {
            "legacy" => false,
            "coset" => true,
            other => return Err(format!("fiber_gen `{other}` is not legacy or coset")),
        };
        let scale = params.u64_or("fiber_scale", 1)?;
        if scale == 0 {
            return Err("fiber_scale must be at least 1".into());
        }
        if scale > 1 && !coset {
            return Err("fiber_scale > 1 needs fiber_gen=coset".into());
        }
        Ok(Self { coset, scale })
    }

    /// The counts a non-default build adds to the set-up record; empty
    /// for the default, so legacy records keep their keys.
    fn counts(self) -> BTreeMap<String, u64> {
        let mut c = BTreeMap::new();
        if self != Self::default() {
            c.insert("fiber_gen_coset".into(), u64::from(self.coset));
            c.insert("fiber_scale".into(), self.scale);
        }
        c
    }
}

/// [`cofactor_fiber`], or the subgroup and closure [`FiberBuild`] names.
pub fn cofactor_fiber_with<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    r: u64,
    h: u64,
    build: FiberBuild,
    ops: &mut GroupOps,
) -> CofactorFiber<G::Elt> {
    if !build.coset {
        return cofactor_fiber(g, fb, r, h, ops);
    }
    let scale = build.scale;
    let mut points = vec![g.identity()];
    let mut keys: HashSet<u64> = HashSet::new();
    keys.insert(g.key(&g.identity()));
    let mut key_probes = 1u64;
    let mut samples = 0u64;
    if h > 1 {
        let mult_by = r.checked_mul(scale).expect("r·d fits a u64");
        for &p in &fb.points {
            if scale == 1 && points.len() as u64 >= h {
                break;
            }
            samples += 1;
            let k = g.mul(ops, p, mult_by);
            key_probes += 1;
            if keys.contains(&g.key(&k)) {
                continue;
            }
            // The cosets `old + jK` for `j` below the index of `old` in
            // `⟨old, K⟩`: the first `jK` already in `old` closes the group.
            let old: Vec<G::Elt> = points.clone();
            let mut mult = k;
            while {
                key_probes += 1;
                !keys.contains(&g.key(&mult))
            } {
                for &o in &old {
                    let q = if g.is_identity(&o) {
                        mult
                    } else {
                        g.add(ops, o, mult)
                    };
                    key_probes += 1;
                    if keys.insert(g.key(&q)) {
                        points.push(q);
                    }
                }
                mult = g.add(ops, mult, k);
                if points.len() as u64 > h {
                    break;
                }
            }
        }
    }
    CofactorFiber {
        complete: if scale == 1 {
            points.len() as u64 == h
        } else {
            samples == fb.points.len() as u64
        },
        points,
        h,
        samples,
        key_probes,
    }
}

const FIBER_GEN_PARAM: (&str, &str) = (
    "fiber_gen",
    "legacy (default; replays committed sessions) or coset (closure stops at the first multiple already in the subgroup)",
);
const FIBER_SCALE_PARAM: (&str, &str) = (
    "fiber_scale",
    "d: use the subgroup [d]E(F_q)[h] instead of the whole fiber (needs fiber_gen=coset; default 1)",
);

/// The lifts `T + K` of a target, the identity's own first.  `|H| − 1`
/// additions, charged.
fn lifts<G: CountedGroup>(
    g: &G,
    fiber: &CofactorFiber<G::Elt>,
    ops: &mut GroupOps,
    t: G::Elt,
) -> Vec<G::Elt> {
    fiber
        .points
        .iter()
        .enumerate()
        .map(|(i, &k)| if i == 0 { t } else { g.add(ops, t, k) })
        .collect()
}

/// Which mode a fiber-aware oracle runs in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiberMode {
    /// The blind oracle: one lift.
    Blind,
    /// Every lift probed or descended.
    Lifts,
    /// A fiber-closed table (table oracles) — one probe.
    Closed,
    /// One combined algebraic system (descent oracle).
    Combined,
}

impl FiberMode {
    pub fn parse(params: &Params, closed_name: &str) -> Result<Self, String> {
        match params.get("fiber").unwrap_or("closed") {
            "blind" => Ok(Self::Blind),
            "lifts" => Ok(Self::Lifts),
            "closed" if closed_name == "closed" => Ok(Self::Closed),
            "combined" if closed_name == "combined" => Ok(Self::Combined),
            other => Err(format!(
                "fiber mode `{other}` is not blind, lifts or {closed_name}"
            )),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Blind => "blind",
            Self::Lifts => "lifts",
            Self::Closed => "closed",
            Self::Combined => "combined",
        }
    }
}

fn fiber_counts<E>(fiber: Option<&CofactorFiber<E>>, build: FiberBuild) -> BTreeMap<String, u64> {
    let mut counts = build.counts();
    if let Some(f) = fiber {
        counts.insert("fiber_size".into(), f.size() as u64);
        counts.insert("fiber_cofactor".into(), f.h);
        counts.insert("fiber_complete".into(), u64::from(f.complete));
        counts.insert("fiber_samples".into(), f.samples);
    }
    counts
}

// ── the fiber-closed pair table ────────────────────────────────────

/// A pair table whose entries are `P_i + P_j + K` for every `K` in the
/// fiber: keyed like [`PairTable`] (the group key without its sign bit),
/// so one probe answers "is `±T` a pair sum up to the fiber".
pub struct FiberPairTable {
    map: HashMap<u64, (u32, u32, u16)>,
    pub kind: &'static str,
    pub entries: u64,
    pub build_ops: GroupOps,
    /// Hash inserts attempted while building, for `bitop/v1` only.
    pub insert_attempts: u64,
}

impl FiberPairTable {
    pub fn build<G: CountedGroup>(
        g: &G,
        fb: &FactorBase<G::Elt>,
        fiber: &CofactorFiber<G::Elt>,
        negation_folded: bool,
    ) -> Self {
        let mut ops = GroupOps::default();
        let mut map: HashMap<u64, (u32, u32, u16)> = HashMap::new();
        let mut attempts = 0u64;
        let mut insert = |i: usize, j: usize, ops: &mut GroupOps| {
            let s = g.add(ops, fb.points[i], fb.points[j]);
            for (k, &kp) in fiber.points.iter().enumerate() {
                let sk = if k == 0 { s } else { g.add(ops, s, kp) };
                if g.is_identity(&sk) {
                    continue;
                }
                attempts += 1;
                map.entry(g.key(&sk) >> 1)
                    .or_insert((i as u32, j as u32, k as u16));
            }
        };
        if negation_folded {
            let groups: Vec<&Vec<usize>> = fb
                .abscissa_list
                .iter()
                .map(|x| fb.x_index.get(x).expect("every listed abscissa is indexed"))
                .collect();
            for ai in 0..groups.len() {
                let i = groups[ai][0];
                for group in groups.iter().skip(ai) {
                    for &j in group.iter() {
                        insert(i, j, &mut ops);
                    }
                }
            }
        } else {
            let n = fb.points.len();
            for i in 0..n {
                for j in i..n {
                    insert(i, j, &mut ops);
                }
            }
        }
        Self {
            kind: if negation_folded {
                "fiber_closed_negation_folded"
            } else {
                "fiber_closed_full"
            },
            entries: map.len() as u64,
            map,
            build_ops: ops,
            insert_attempts: attempts,
        }
    }

    /// Indices `(i, j)` with `P_i + P_j ∈ ±target − H`, if the table holds
    /// such a pair: one probe, at most two additions.
    pub fn probe<G: CountedGroup>(
        &self,
        g: &G,
        fb: &FactorBase<G::Elt>,
        fiber: &CofactorFiber<G::Elt>,
        ops: &mut GroupOps,
        lookups: &mut u64,
        target: G::Elt,
    ) -> Option<(usize, usize)> {
        *lookups += 1;
        let &(i, j, k) = self.map.get(&(g.key(&target) >> 1))?;
        let (i, j, k) = (i as usize, j as usize, k as usize);
        let mut s = g.add(ops, fb.points[i], fb.points[j]);
        if k != 0 {
            s = g.add(ops, s, fiber.points[k]);
        }
        if s == target {
            Some((i, j))
        } else if s == g.neg(target) {
            Some((fb.neg_index[i], fb.neg_index[j]))
        } else {
            None
        }
    }
}

// ── the generic fiber-aware MITM oracle ────────────────────────────

/// Meet in the middle, fiber-aware: `mitm-fiber:fiber=closed|lifts|blind`,
/// `negation_folded=1` as for `mitm`.
pub struct FiberMitmOracle<E> {
    m: u32,
    mode: FiberMode,
    fiber: Option<CofactorFiber<E>>,
    build: FiberBuild,
    plain: Option<PairTable>,
    closed: Option<FiberPairTable>,
}

impl<E> FiberMitmOracle<E> {
    pub fn new(m: u32) -> Self {
        Self {
            m,
            mode: FiberMode::Closed,
            fiber: None,
            build: FiberBuild::default(),
            plain: None,
            closed: None,
        }
    }
}

impl<G: CountedGroup> DecompositionOracle<G> for FiberMitmOracle<G::Elt> {
    fn name(&self) -> &str {
        "mitm-fiber"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "meet in the middle over the cofactor fiber ({}), {}",
            params.get("fiber").unwrap_or("closed"),
            if params.flag("negation_folded") {
                "table built with the negation fold"
            } else {
                "table built in full"
            }
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("fiber", "closed (table holds P_i + P_j + K, one probe), lifts (h probes per target) or blind"),
            ("negation_folded", "1 to build one of each ±pair"),
            FIBER_GEN_PARAM,
            FIBER_SCALE_PARAM,
        ]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        self.mode = FiberMode::parse(params, "closed")?;
        self.build = FiberBuild::parse(params)?;
        let fiber = cofactor_fiber_with(ctx.group, fb, ctx.r, ctx.cofactor, self.build, ops);
        let folded = params.flag("negation_folded");
        match self.mode {
            FiberMode::Closed => {
                let t = FiberPairTable::build(ctx.group, fb, &fiber, folded);
                ops.merge(t.build_ops);
                self.closed = Some(t);
            }
            _ => {
                let t = if folded {
                    PairTable::build_negation_folded(ctx.group, fb)
                } else {
                    PairTable::build(ctx.group, fb)
                };
                ops.merge(t.build_ops);
                self.plain = Some(t);
            }
        }
        self.fiber = Some(fiber);
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut c = fiber_counts(self.fiber.as_ref(), self.build);
        if let Some(t) = &self.closed {
            c.insert("fiber_table_entries".into(), t.entries);
        }
        if let Some(t) = &self.plain {
            c.insert("pair_table_entries".into(), t.entries);
        }
        c
    }

    fn setup_work_native(&self) -> BTreeMap<String, u64> {
        let mut c = match (&self.closed, &self.plain) {
            (Some(t), _) => {
                BTreeMap::from([("table_insert_attempts".to_string(), t.insert_attempts)])
            }
            (None, Some(t)) => t.work_counters(),
            (None, None) => BTreeMap::new(),
        };
        add_fiber_key_probes(&mut c, self.fiber.as_ref());
        c.insert("counted:table_insert_attempts".to_string(), 1);
        c
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let g = ctx.group;
        let fiber = self.fiber.as_ref()?;
        match self.mode {
            FiberMode::Closed => {
                let table = self.closed.as_ref()?;
                match self.m {
                    2 => table
                        .probe(g, fb, fiber, ops, &mut counters.lookups, point)
                        .map(|(i, j)| vec![i, j]),
                    3 => {
                        for i in 0..fb.points.len() {
                            let s = g.add(ops, point, g.neg(fb.points[i]));
                            if g.is_identity(&s) {
                                continue;
                            }
                            if let Some((j, k)) =
                                table.probe(g, fb, fiber, ops, &mut counters.lookups, s)
                            {
                                return Some(vec![i, j, k]);
                            }
                        }
                        None
                    }
                    _ => None,
                }
            }
            FiberMode::Blind | FiberMode::Lifts => {
                let table = self.plain.as_ref()?;
                let targets = if self.mode == FiberMode::Blind {
                    vec![point]
                } else {
                    lifts(g, fiber, ops, point)
                };
                for t in targets {
                    let hit = crate::cryptanalysis::ic_boundary::decompose_mitm(
                        g, fb, table, self.m, ops, counters, t,
                    );
                    if hit.is_some() {
                        return hit;
                    }
                }
                None
            }
            FiberMode::Combined => None,
        }
    }
}

// ── the Frobenius-folded fiber table (Koblitz) ─────────────────────

/// Bits of a table entry that hold the second index.
const FROB_INDEX_BITS: u32 = 26;

/// [`FrobeniusPairTable`] with the fiber folded in: entries
/// `σ^t(P_rep + P_j) + K` for every `K ∈ E(F_q)[h]`.  Needs `H ⊂ E(F_2)`
/// (so that the Frobenius fixes every `K` and the rotation argument of
/// the plain table carries over), which holds on every Koblitz instance
/// and is checked in `build`.
pub struct FiberFrobeniusPairTable {
    map: HashMap<u64, (u32, u32, u16)>,
    canon: FrobeniusCanon,
    member_of: Vec<u32>,
    n: u32,
    pub entries: u64,
    pub build_ops: GroupOps,
    pub build_canonicalisations: u64,
    /// Hash inserts attempted while building, for `bitop/v1` only.
    pub insert_attempts: u64,
}

impl FiberFrobeniusPairTable {
    pub fn build(
        inst: &BinaryInstance,
        fb: &FactorBase<FastPoint>,
        fiber: &CofactorFiber<FastPoint>,
    ) -> Result<Self, String> {
        let g = BinaryGroup(&inst.fast);
        let n = inst.n as usize;
        let f = fb.points.len();
        if f >= (1usize << FROB_INDEX_BITS) {
            return Err("base too large for the folded table's index".into());
        }
        for &k in &fiber.points {
            if inst.fast.frobenius_k(k, 1) != k {
                return Err("a fiber point is not fixed by the Frobenius: the folded fiber table needs H ⊂ E(F_2)".into());
            }
        }
        let canon = FrobeniusCanon::new(&inst.gf, inst.n).ok_or("no normal basis found")?;
        let mut member_of = vec![0u32; f * n];
        for i in 0..f {
            let mut p = fb.points[i];
            member_of[i * n] = i as u32;
            for t in 1..n {
                p = inst.fast.frobenius_k(p, 1);
                let idx = *fb
                    .point_index
                    .get(&p.pack())
                    .ok_or("the base is not closed under the Frobenius")?;
                member_of[i * n + t] = idx as u32;
            }
        }
        let mut rep_of_col = vec![usize::MAX; fb.columns];
        for i in 0..f {
            let c = fb.col_of[i];
            if rep_of_col[c] == usize::MAX {
                rep_of_col[c] = i;
            }
        }
        let mut ops = GroupOps::default();
        let mut canons = 0u64;
        let mut attempts = 0u64;
        let mut map: HashMap<u64, (u32, u32, u16)> = HashMap::new();
        for (col, &rep) in rep_of_col.iter().enumerate() {
            if rep == usize::MAX {
                continue;
            }
            for j in 0..f {
                if fb.col_of[j] < col {
                    continue;
                }
                let u = g.add(&mut ops, fb.points[rep], fb.points[j]);
                for (k, &kp) in fiber.points.iter().enumerate() {
                    let uk = if k == 0 { u } else { g.add(&mut ops, u, kp) };
                    if g.is_identity(&uk) {
                        continue;
                    }
                    canons += 1;
                    attempts += 1;
                    let (key, t_u) = canon.canon_with_shift(uk.x);
                    map.entry(key).or_insert((
                        rep as u32,
                        (j as u32) | (t_u << FROB_INDEX_BITS),
                        k as u16,
                    ));
                }
            }
        }
        Ok(Self {
            entries: map.len() as u64,
            map,
            canon,
            member_of,
            n: inst.n,
            build_ops: ops,
            build_canonicalisations: canons,
            insert_attempts: attempts,
        })
    }

    pub fn probe(
        &self,
        g: &BinaryGroup,
        fb: &FactorBase<FastPoint>,
        fiber: &CofactorFiber<FastPoint>,
        ops: &mut GroupOps,
        ctr: &mut OracleCounters,
        target: FastPoint,
    ) -> Option<(usize, usize)> {
        ctr.canonicalisations += 1;
        let (key, t_r) = self.canon.canon_with_shift(target.x);
        ctr.lookups += 1;
        let &(i, packed, k) = self.map.get(&key)?;
        let n = self.n;
        let j = (packed & ((1u32 << FROB_INDEX_BITS) - 1)) as usize;
        let t_u = packed >> FROB_INDEX_BITS;
        let t = ((t_u + n - t_r) % n) as usize;
        let i2 = self.member_of[i as usize * n as usize + t] as usize;
        let j2 = self.member_of[j * n as usize + t] as usize;
        let mut s = g.add(ops, fb.points[i2], fb.points[j2]);
        if k != 0 {
            s = g.add(ops, s, fiber.points[k as usize]);
        }
        if s == target {
            Some((i2, j2))
        } else if s == g.neg(target) {
            Some((fb.neg_index[i2], fb.neg_index[j2]))
        } else {
            ctr.frobfold_mismatches += 1;
            None
        }
    }
}

/// `mitm-frobenius-fiber:fiber=closed|lifts|blind` on a Koblitz instance.
pub struct FiberFrobeniusMitmOracle<'i> {
    m: u32,
    instance: &'i BinaryInstance,
    mode: FiberMode,
    fiber: Option<CofactorFiber<FastPoint>>,
    build: FiberBuild,
    plain: Option<FrobeniusPairTable>,
    closed: Option<FiberFrobeniusPairTable>,
}

impl<'i> FiberFrobeniusMitmOracle<'i> {
    pub fn new(m: u32, instance: &'i BinaryInstance) -> Self {
        Self {
            m,
            instance,
            mode: FiberMode::Closed,
            fiber: None,
            build: FiberBuild::default(),
            plain: None,
            closed: None,
        }
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for FiberFrobeniusMitmOracle<'_> {
    fn name(&self) -> &str {
        "mitm-frobenius-fiber"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "a pair table folded by the signed Frobenius, fiber-aware ({})",
            params.get("fiber").unwrap_or("closed")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("fiber", "closed, lifts or blind"),
            FIBER_GEN_PARAM,
            FIBER_SCALE_PARAM,
        ]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        self.mode = FiberMode::parse(params, "closed")?;
        self.build = FiberBuild::parse(params)?;
        let fiber = cofactor_fiber_with(ctx.group, fb, ctx.r, ctx.cofactor, self.build, ops);
        match self.mode {
            FiberMode::Closed => {
                let t = FiberFrobeniusPairTable::build(self.instance, fb, &fiber)?;
                ops.merge(t.build_ops);
                self.closed = Some(t);
            }
            _ => {
                let t = FrobeniusPairTable::build(self.instance, fb)
                    .ok_or("this instance has no Koblitz structure to fold by")?;
                ops.merge(t.build_ops);
                self.plain = Some(t);
            }
        }
        self.fiber = Some(fiber);
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        let mut c = fiber_counts(self.fiber.as_ref(), self.build);
        if let Some(t) = &self.closed {
            c.insert("fiber_table_entries".into(), t.entries);
            c.insert("canonicalisations".into(), t.build_canonicalisations);
        }
        if let Some(t) = &self.plain {
            c.insert("pair_table_entries".into(), t.entries);
        }
        c
    }

    fn setup_work_native(&self) -> BTreeMap<String, u64> {
        let mut c = match (&self.closed, &self.plain) {
            (Some(t), _) => {
                // The table builds its own canonicalisation tables.
                let canon = t.canon.build_work();
                let mut c = BTreeMap::new();
                for (name, v) in [
                    ("table_insert_attempts", t.insert_attempts),
                    ("canon_build_word_ops", canon.word_ops),
                    ("canon_build_field_sqrs", canon.field_sqrs),
                ] {
                    if v > 0 {
                        c.insert(name.to_string(), v);
                    }
                }
                c
            }
            (None, Some(t)) => t.work_counters(),
            (None, None) => BTreeMap::new(),
        };
        add_fiber_key_probes(&mut c, self.fiber.as_ref());
        c.insert("counted:table_insert_attempts".to_string(), 1);
        c
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: FastPoint,
    ) -> Option<Vec<usize>> {
        let g = ctx.group;
        let fiber = self.fiber.as_ref()?;
        match self.mode {
            FiberMode::Closed => {
                let table = self.closed.as_ref()?;
                match self.m {
                    2 => table
                        .probe(g, fb, fiber, ops, counters, point)
                        .map(|(i, j)| vec![i, j]),
                    3 => {
                        for i in 0..fb.points.len() {
                            let s = g.add(ops, point, g.neg(fb.points[i]));
                            if g.is_identity(&s) {
                                continue;
                            }
                            if let Some((j, k)) = table.probe(g, fb, fiber, ops, counters, s) {
                                return Some(vec![i, j, k]);
                            }
                        }
                        None
                    }
                    _ => None,
                }
            }
            FiberMode::Blind | FiberMode::Lifts => {
                let table = self.plain.as_ref()?;
                let targets = if self.mode == FiberMode::Blind {
                    vec![point]
                } else {
                    lifts(g, fiber, ops, point)
                };
                for t in targets {
                    let hit = crate::cryptanalysis::ic_boundary::decompose_mitm_frobenius(
                        g, fb, table, self.m, ops, counters, t,
                    );
                    if hit.is_some() {
                        return hit;
                    }
                }
                None
            }
            FiberMode::Combined => None,
        }
    }
}

// ── the fiber-aware algebraic oracle ───────────────────────────────

/// `descent-algebraic-fiber:fiber=combined|lifts|blind` with a `solver`,
/// on a subspace base.  `lifts` descends every lift (`h` systems per
/// target); `combined` descends one system with the target abscissa free
/// and the fiber polynomial adjoined ([`descend_fiber`]).
pub struct FiberDescentOracle<'i> {
    summands: u32,
    instance: &'i BinaryInstance,
    v_basis: Vec<u64>,
    harness: SolverHarness,
    totals: SolverTotals,
    mode: FiberMode,
    fiber: Option<CofactorFiber<FastPoint>>,
    build: FiberBuild,
    pub unliftable: u64,
}

impl<'i> FiberDescentOracle<'i> {
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
            mode: FiberMode::Combined,
            fiber: None,
            build: FiberBuild::default(),
            unliftable: 0,
        }
    }

    fn descend_plain(&self, x_r: u64) -> SymbolicDescent {
        descend(
            &self.instance.gf,
            self.instance.b,
            x_r,
            &self.v_basis,
            self.summands,
        )
        .expect("prepare checked the dimension")
    }

    fn descend_combined(&self, fiber_xs: &[u64]) -> SymbolicDescent {
        descend_fiber(
            &self.instance.gf,
            self.instance.b,
            fiber_xs,
            &self.v_basis,
            self.summands,
        )
        .expect("prepare checked the dimension")
    }

    /// Solve one system and lift its solutions against the given lifts.
    fn solve_and_lift(
        &mut self,
        g: &BinaryGroup,
        fb: &FactorBase<FastPoint>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        sys: &SymbolicDescent,
        targets_by_x: &dyn Fn(u64) -> Vec<FastPoint>,
        all_targets: &[FastPoint],
    ) -> Option<Vec<usize>> {
        let system = BooleanSystem {
            equations: sys.equations.clone(),
            n_vars: sys.n_vars,
        };
        let shape = system.shape();
        let verdict = self.harness.solve(&system);
        let exceeded = matches!(verdict, SolverVerdict::BudgetExceeded);
        self.totals.absorb(&shape, self.harness.cost(), exceeded);
        // The descent (counted on this thread since `decompose` began),
        // the copy, and two shape computations (here and in the harness).
        charge_system_build(&mut self.totals.extra, &system, &shape, 2);
        let SolverVerdict::Solved(solutions) = verdict else {
            return None;
        };
        for v in solutions {
            let xs = sys.lift(v);
            // The combined system names the lift it found; the plain one
            // is checked against every lift it was asked about.
            let candidates: Vec<FastPoint> = if self.mode == FiberMode::Combined {
                // `lift_target` reads `n` bits (a word op each); the lookup
                // compares the abscissa of every lift (a field add each).
                *self
                    .totals
                    .extra
                    .entry("bitop_word_ops".to_string())
                    .or_insert(0) += u64::from(sys.n);
                *self
                    .totals
                    .extra
                    .entry("bitop_field_adds".to_string())
                    .or_insert(0) += all_targets.len() as u64;
                targets_by_x(sys.lift_target(v))
            } else {
                all_targets.to_vec()
            };
            let mut lift_pending = true;
            for t in candidates {
                let lifted = lift_abscissae(g, fb, ops, &xs, t);
                // The lift's share of the descent counter is charged once,
                // with the first candidate; each candidate probes the base.
                if lift_pending {
                    charge_lift(&mut self.totals.extra, fb, &xs, lifted.is_some());
                    lift_pending = false;
                } else {
                    *self
                        .totals
                        .extra
                        .entry("bitop_hash_probes".to_string())
                        .or_insert(0) += lift_probes(fb, &xs, lifted.is_some());
                }
                if let Some(indices) = lifted {
                    return Some(indices);
                }
            }
            if lift_pending {
                // No candidate: the lift itself is still work done.
                crate::cryptanalysis::pq_descent_symbolic::take_descent_work()
                    .charge(&mut self.totals.extra);
            }
            counters.lift_failures += 1;
        }
        self.unliftable += 1;
        counters.unliftable_systems += 1;
        None
    }
}

impl<'a> DecompositionOracle<BinaryGroup<'a>> for FiberDescentOracle<'_> {
    fn name(&self) -> &str {
        "descent-algebraic-fiber"
    }

    fn summands(&self) -> u32 {
        self.summands
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "Weil-descend S_{} over the base's abscissa subspace, fiber-aware ({}), solved with {}",
            self.summands + 1,
            params.get("fiber").unwrap_or("combined"),
            self.harness.solver.name()
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("m", "summands, 2 (descends S3) or 3 (descends S4)"),
            ("fiber", "combined (one system, target abscissa free, fiber polynomial adjoined), lifts (one system per lift) or blind"),
            FIBER_GEN_PARAM,
            FIBER_SCALE_PARAM,
        ]
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<BinaryGroup<'a>>,
        fb: &FactorBase<FastPoint>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        self.mode = FiberMode::parse(params, "combined")?;
        if !(2..=3).contains(&self.summands) {
            return Err(format!(
                "descent-algebraic-fiber takes m = 2 or 3, got {}",
                self.summands
            ));
        }
        let basis = fb.subspace_basis.as_ref().ok_or(
            "descent-algebraic-fiber needs a factor base whose abscissae form an F_2-subspace",
        )?;
        let cap = max_n_prime(self.summands);
        if basis.len() as u32 > cap {
            return Err(format!(
                "the descent's monomial mask holds n' ≤ {cap} at m = {}; this base has n' = {}",
                self.summands,
                basis.len()
            ));
        }
        if self.mode == FiberMode::Combined
            && self.summands * basis.len() as u32 + self.instance.n > 64
        {
            return Err(format!(
                "the combined descent needs m·n' + n = {} boolean variables; the mask holds 64",
                self.summands * basis.len() as u32 + self.instance.n
            ));
        }
        self.v_basis = basis.clone();
        self.build = FiberBuild::parse(params)?;
        let fiber = cofactor_fiber_with(ctx.group, fb, ctx.r, ctx.cofactor, self.build, ops);
        let probe_xs: Vec<u64> = fiber
            .points
            .iter()
            .map(|k| ctx.group.add(ops, ctx.generator, *k).x)
            .collect();
        let _ = crate::cryptanalysis::pq_descent_symbolic::take_descent_work();
        let probe = match self.mode {
            FiberMode::Combined => BooleanSystem {
                equations: self.descend_combined(&probe_xs).equations,
                n_vars: self.summands as usize * basis.len() + self.instance.n as usize,
            },
            _ => BooleanSystem {
                equations: self.descend_plain(self.instance.generator.x).equations,
                n_vars: self.summands as usize * basis.len(),
            },
        };
        let shape = probe.shape();
        // The probe system moves (no copy): its descent and one shape.
        crate::cryptanalysis::pq_descent_symbolic::take_descent_work()
            .charge(&mut self.totals.extra);
        *self
            .totals
            .extra
            .entry("bitop_word_ops".to_string())
            .or_insert(0) += probe.shape_word_ops(&shape);
        if !self.harness.solver.accepts(&shape) {
            return Err(format!(
                "solver `{}` declines a system of {} unknowns and {} equations",
                self.harness.solver.name(),
                shape.n_vars,
                shape.n_equations
            ));
        }
        self.fiber = Some(fiber);
        Ok(())
    }

    fn setup_native(&self) -> BTreeMap<String, u64> {
        fiber_counts(self.fiber.as_ref(), self.build)
    }

    fn setup_work_native(&self) -> BTreeMap<String, u64> {
        // No pair table; the fiber's generation probes its key set.
        let mut c = BTreeMap::new();
        add_fiber_key_probes(&mut c, self.fiber.as_ref());
        c.insert("counted:table_insert_attempts".to_string(), 1);
        c
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
        let g = ctx.group;
        let fiber = self.fiber.clone()?;
        // Work left on this thread's descent counter by anything else is
        // not this target's.
        let _ = crate::cryptanalysis::pq_descent_symbolic::take_descent_work();
        let all: Vec<FastPoint> = match self.mode {
            FiberMode::Blind => vec![point],
            _ => lifts(g, &fiber, ops, point),
        };
        match self.mode {
            FiberMode::Combined => {
                let xs: Vec<u64> = all.iter().map(|t| t.x).collect();
                let sys = self.descend_combined(&xs);
                let by_x = |x: u64| -> Vec<FastPoint> {
                    all.iter().copied().filter(|t| t.x == x).collect()
                };
                self.solve_and_lift(g, fb, ops, counters, &sys, &by_x, &all)
            }
            _ => {
                for t in &all {
                    let sys = self.descend_plain(t.x);
                    let one = [*t];
                    let by_x = |_x: u64| -> Vec<FastPoint> { vec![*t] };
                    if let Some(hit) = self.solve_and_lift(g, fb, ops, counters, &sys, &by_x, &one)
                    {
                        return Some(hit);
                    }
                }
                None
            }
        }
    }

    fn last_system(&self) -> Option<SystemShape> {
        self.totals.shape.clone()
    }

    fn solver_totals(&self) -> Option<SolverTotals> {
        Some(self.totals.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::{binary_subspace_factor_base, koblitz_instance};

    #[test]
    fn the_cofactor_fiber_of_a_koblitz_curve_is_its_rational_torsion() {
        let inst = koblitz_instance(0, 13).expect("the a = 0 Koblitz curve over F_{2^13}");
        let basis: Vec<u64> = (0..5).map(|k| 1u64 << k).collect();
        let fb = binary_subspace_factor_base(&inst, &basis);
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let fiber = cofactor_fiber(&g, &fb, inst.r, inst.cofactor, &mut ops);
        assert_eq!(fiber.h, inst.cofactor);
        assert!(fiber.complete, "found {} of {}", fiber.size(), fiber.h);
        for &k in &fiber.points {
            assert!(g.mul(&mut ops, k, inst.cofactor).infinity);
            assert_eq!(inst.fast.frobenius_k(k, 1), k, "the h-torsion is over F_2");
        }
    }

    #[test]
    fn the_fiber_closed_table_finds_every_lift_the_plain_table_finds() {
        let inst = koblitz_instance(0, 13).expect("the a = 0 Koblitz curve over F_{2^13}");
        let basis: Vec<u64> = (0..5).map(|k| 1u64 << k).collect();
        let fb = binary_subspace_factor_base(&inst, &basis);
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let fiber = cofactor_fiber(&g, &fb, inst.r, inst.cofactor, &mut ops);
        let plain = PairTable::build_negation_folded(&g, &fb);
        let closed = FiberPairTable::build(&g, &fb, &fiber, true);
        assert!(closed.entries > plain.entries);
        let mut lookups = 0u64;
        let mut hits_plain = 0;
        let mut hits_closed = 0;
        let mut p = inst.generator;
        for _ in 0..2000 {
            p = g.add(&mut ops, p, inst.generator);
            let mut any_plain = false;
            for &k in &fiber.points {
                let t = g.add(&mut ops, p, k);
                if plain.probe(&g, &fb, &mut ops, &mut lookups, t).is_some() {
                    any_plain = true;
                }
            }
            let c = closed.probe(&g, &fb, &fiber, &mut ops, &mut lookups, p);
            if let Some((i, j)) = c {
                let s = g.add(&mut ops, fb.points[i], fb.points[j]);
                let d = g.add(&mut ops, s, g.neg(p));
                assert!(
                    g.mul(&mut ops, d, inst.cofactor).infinity,
                    "a closed hit is a lift"
                );
            }
            hits_plain += usize::from(any_plain);
            hits_closed += usize::from(c.is_some());
        }
        assert!(
            hits_closed >= hits_plain,
            "closed {hits_closed} < plain-over-lifts {hits_plain}"
        );
        assert!(hits_closed > 0);
    }

    #[test]
    fn the_combined_descent_names_a_lift() {
        use crate::cryptanalysis::ic_framework::solvers::solver_by_name;
        let inst = koblitz_instance(0, 13).expect("the a = 0 Koblitz curve over F_{2^13}");
        let basis: Vec<u64> = (0..5).map(|k| 1u64 << k).collect();
        let fb = binary_subspace_factor_base(&inst, &basis);
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let fiber = cofactor_fiber(&g, &fb, inst.r, inst.cofactor, &mut ops);
        // A decomposable lift: P_1 + P_2 + K for base points and a fiber point.
        let p1 = fb.points[3];
        let p2 = fb.points[7];
        let k = fiber.points[1];
        let p12 = g.add(&mut ops, p1, p2);
        let s = g.add(&mut ops, p12, g.neg(k));
        // The target is the fiber element with K = O: s itself; its fiber holds s + K = P_1 + P_2.
        let all = lifts(&g, &fiber, &mut ops, s);
        let xs: Vec<u64> = all.iter().map(|t| t.x).collect();
        let sys = descend_fiber(&inst.gf, inst.b, &xs, &basis, 2).unwrap();
        assert_eq!(sys.n_vars, 2 * 5 + 13);
        let solver = solver_by_name("f4-f2").unwrap();
        let (verdict, _) = solver.solve(
            &BooleanSystem {
                equations: sys.equations.clone(),
                n_vars: sys.n_vars,
            },
            &Params::default(),
            None,
        );
        let SolverVerdict::Solved(sols) = verdict else {
            panic!("no solution")
        };
        let mut found = false;
        for v in sols {
            let abs = sys.lift(v);
            let x = sys.lift_target(v);
            let targets: Vec<FastPoint> = all.iter().copied().filter(|t| t.x == x).collect();
            assert!(
                !targets.is_empty(),
                "the named target abscissa is in the fiber"
            );
            for t in targets {
                if lift_abscissae(&g, &fb, &mut ops, &abs, t).is_some() {
                    found = true;
                }
            }
        }
        assert!(found, "the planted lift was recovered");
    }

    #[test]
    fn the_coset_closure_builds_the_same_fiber_for_fewer_additions() {
        let inst = koblitz_instance(0, 13).expect("the a = 0 Koblitz curve over F_{2^13}");
        let basis: Vec<u64> = (0..5).map(|k| 1u64 << k).collect();
        let fb = binary_subspace_factor_base(&inst, &basis);
        let g = BinaryGroup(&inst.fast);
        let mut legacy_ops = GroupOps::default();
        let legacy = cofactor_fiber(&g, &fb, inst.r, inst.cofactor, &mut legacy_ops);
        let mut coset_ops = GroupOps::default();
        let build = FiberBuild {
            coset: true,
            scale: 1,
        };
        let coset = cofactor_fiber_with(&g, &fb, inst.r, inst.cofactor, build, &mut coset_ops);
        assert!(coset.complete);
        let key_set = |f: &CofactorFiber<FastPoint>| -> HashSet<u64> {
            f.points.iter().map(|k| g.key(k)).collect()
        };
        assert_eq!(key_set(&legacy), key_set(&coset));
        assert!(coset_ops.gae() <= legacy_ops.gae());
        let mut legacy_again = GroupOps::default();
        let again = cofactor_fiber_with(
            &g,
            &fb,
            inst.r,
            inst.cofactor,
            FiberBuild::default(),
            &mut legacy_again,
        );
        assert_eq!(
            legacy_again.gae(),
            legacy_ops.gae(),
            "the default is the legacy count"
        );
        assert_eq!(again.samples, legacy.samples);
    }

    #[test]
    fn a_scaled_fiber_is_the_image_of_multiplication_by_d() {
        let inst = koblitz_instance(0, 13).expect("the a = 0 Koblitz curve over F_{2^13}");
        let basis: Vec<u64> = (0..5).map(|k| 1u64 << k).collect();
        let fb = binary_subspace_factor_base(&inst, &basis);
        let g = BinaryGroup(&inst.fast);
        let mut ops = GroupOps::default();
        let full = cofactor_fiber(&g, &fb, inst.r, inst.cofactor, &mut ops);
        for d in [1u64, 2, 4] {
            let build = FiberBuild {
                coset: true,
                scale: d,
            };
            let sub = cofactor_fiber_with(&g, &fb, inst.r, inst.cofactor, build, &mut ops);
            let mut image: HashSet<u64> = HashSet::new();
            for &k in &full.points {
                image.insert(g.key(&g.mul(&mut ops, k, d)));
            }
            let got: HashSet<u64> = sub.points.iter().map(|k| g.key(k)).collect();
            assert_eq!(got, image, "d = {d}");
        }
    }

    #[test]
    fn fiber_build_parameters_are_checked() {
        let p = |kv: &[(&str, &str)]| {
            let mut ps = Params::default();
            for (k, v) in kv {
                ps.set(k, *v);
            }
            FiberBuild::parse(&ps)
        };
        assert_eq!(p(&[]).unwrap(), FiberBuild::default());
        assert!(p(&[("fiber_gen", "coset")]).unwrap().coset);
        assert!(p(&[("fiber_scale", "2")]).is_err());
        assert!(p(&[("fiber_gen", "coset"), ("fiber_scale", "0")]).is_err());
        assert!(p(&[("fiber_gen", "other")]).is_err());
        assert!(FiberBuild::default().counts().is_empty());
    }
}
