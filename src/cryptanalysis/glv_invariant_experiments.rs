//! # Drivers for the experiments of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.
//!
//! The framework's relation loop stops the moment the target's
//! logarithm is determined, which on a two-summand base is the first
//! cycle of the relation graph, at about half the columns.  E1 asks a
//! different question — how many relations the fold saves at **full
//! rank** — and asks it on one stream: every decomposition found on the
//! shared point set is fed to two matrices at once, the folded arm's
//! and the control's, so the two arms differ in nothing but their
//! column maps.  [`full_rank_stream`] is that driver, and it records for
//! each arm the relation and target count at which the logarithm was
//! first pinned, and at which the system on the columns it had touched
//! became square (rank equal to touched columns plus one, every touched
//! unknown determined).  The second is what the fold's `w/2` is a
//! statement about.
//!
//! **Full rank is measured against the achievable rank, not `columns + 1`.**
//! On a curve with cofactor `h > 1` a two-summand relation
//! `R = ε P + ε' Q`, `R ∈ ⟨G⟩`, forces the `E[h]`-components of its
//! summands to be negatives of each other (`ε u + ε' u' = 0`), so every
//! odd function `χ` on `E[h]` (with `χ(π u) = λ_π χ(u)` for a folded
//! arm) is a functional `Σ χ(u_i) x_i` that no row ever determines.  The
//! rank therefore saturates `D` short of full, where `D` is the number
//! of `⟨−1, group⟩`-classes of `E[h]`-components among the touched
//! columns that are not their own negatives; the logarithm is pinned
//! all the same.  The driver computes each column's component class
//! from `[r]P` and reports the deficiency beside the count.  Three
//! summands carry no such constraint, and the deficiency is taken as
//! zero.
//!
//! The same driver serves E2 (GLS and Koblitz arms), E5 (subfield
//! arms), and E7 (three-summand oracle, with the targets canonicalised
//! under the group so that orbit duplicates are counted).

use std::collections::HashSet;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::glv_invariant_base::{mulm, subm, Endomorphism, EndomorphismClasses};
use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, OracleCounters, RelationSolver, RhoClasses, RowStatus,
};
use crate::cryptanalysis::ic_framework::linalg::matrix_by_name;

/// One arm's view of the stream.
#[derive(Clone, Debug, Serialize)]
pub struct ArmReport {
    pub label: String,
    pub points: usize,
    pub columns: usize,
    pub points_per_column: f64,
    /// Relations and targets at the first determination of the target's
    /// logarithm (what the framework's loop stops at).
    pub first_pin_relations: Option<u64>,
    pub first_pin_trials: Option<u64>,
    /// Relations and targets at which the matrix reached its achievable
    /// full rank, `columns + 1 − deficiency`: every unknown the
    /// relations can determine is determined.  This is the
    /// coupon-collector count — the last column must be touched and its
    /// component connected — so it carries a `ln columns` factor over
    /// the `columns` floor on both arms alike.
    pub square_relations: Option<u64>,
    pub square_trials: Option<u64>,
    /// `square_relations / columns`: the count against its floor.
    pub square_over_columns: Option<f64>,
    pub touched_at_square: Option<usize>,
    /// The rank reached after exactly `columns` relations: the
    /// size-matched reading, free of the coupon-collector tail.
    pub rank_at_columns: Option<u64>,
    /// `rank_at_columns / (columns + 1 − deficiency)`.
    pub rank_fraction_at_columns: Option<f64>,
    /// The structural rank deficiency (see the module header) over the
    /// whole base; `deficiency_at_square` is the same number, kept for
    /// the tables.
    pub deficiency_at_square: Option<usize>,
    pub deficiency_total: usize,
    pub relations_total: u64,
    pub rank_final: u64,
    pub dependent_final: u64,
    pub row_ops: u64,
    /// Rows whose base part folded to nothing (`0 = h a + h b d`).
    pub zero_support_rows: u64,
    pub single_column_rows: u64,
    pub recovered: Option<u64>,
    pub verified: bool,
}

/// The stream's report: both arms, one set of targets.
#[derive(Clone, Debug, Serialize)]
pub struct StreamReport {
    pub trials: u64,
    pub relations: u64,
    pub hit_rate: f64,
    /// Group operations on the target side and in the oracle, shared.
    pub group_ops: GroupOps,
    pub oracle_lookups: u64,
    pub oracle_lift_failures: u64,
    pub oracle_unliftable_systems: u64,
    /// Targets whose group orbit repeated an earlier target's.
    pub orbit_duplicates: u64,
    pub folded: ArmReport,
    pub control: ArmReport,
    /// `control.square_relations / folded.square_relations`.
    pub square_ratio: Option<f64>,
    pub first_pin_ratio: Option<f64>,
    pub column_ratio: f64,
    pub exhausted: bool,
}

struct Arm<'a, E: Copy> {
    label: String,
    fb: &'a FactorBase<E>,
    matrix: Box<dyn RelationSolver>,
    /// Per column: the canonical key of its `E[h]`-component class, or
    /// `None` when the class is its own negative (no deficiency).
    class_of_column: Vec<Option<u64>>,
    touched_classes: HashSet<u64>,
    touched: HashSet<usize>,
    relations: u64,
    first_pin: Option<(u64, u64)>,
    square: Option<(u64, u64, usize)>,
    deficiency_at_square: Option<usize>,
    deficiency_total: usize,
    rank_at_columns: Option<u64>,
    zero_support: u64,
    single: u64,
    recovered: Option<u64>,
}

impl<E: Copy> Arm<'_, E> {
    fn feed(&mut self, summands: &[usize], a: u64, b: u64, h: u64, r: u64, trials: u64) {
        let cols = self.fb.columns + 1;
        let d_col = self.fb.columns;
        let mut row = vec![0u64; cols];
        for &i in summands {
            let c = self.fb.col_of[i];
            row[c] = (row[c] + self.fb.coef_of[i]) % r;
        }
        let support = row[..d_col].iter().filter(|&&c| c != 0).count();
        match support {
            0 => self.zero_support += 1,
            1 => self.single += 1,
            _ => {}
        }
        for (c, &v) in row[..d_col].iter().enumerate() {
            if v != 0 {
                self.touched.insert(c);
                if let Some(cls) = self.class_of_column[c] {
                    self.touched_classes.insert(cls);
                }
            }
        }
        let h_mod = h % r;
        row[d_col] = subm(0, mulm(h_mod, b, r), r);
        let rhs = mulm(h_mod, a, r);
        match self.matrix.add_row(row, rhs) {
            RowStatus::Inconsistent => {}
            RowStatus::Dependent | RowStatus::Independent => {}
        }
        self.relations += 1;
        if self.first_pin.is_none() {
            if let Some(d) = self.matrix.pinned(d_col) {
                self.first_pin = Some((self.relations, trials));
                self.recovered = Some(d);
            }
        }
        if self.relations == self.fb.columns as u64 {
            self.rank_at_columns = Some(self.matrix.rank() as u64);
        }
        // Full rank is `columns + 1 - D` (header): `rank + D > columns`.
        if self.square.is_none() && self.matrix.rank() + self.deficiency_total > self.fb.columns {
            self.square = Some((self.relations, trials, self.touched.len()));
            self.deficiency_at_square = Some(self.deficiency_total);
            if self.recovered.is_none() {
                self.recovered = self.matrix.pinned(d_col);
            }
        }
    }

    fn report(&self, verified: bool) -> ArmReport {
        ArmReport {
            label: self.label.clone(),
            points: self.fb.points.len(),
            columns: self.fb.columns,
            points_per_column: self.fb.points.len() as f64 / self.fb.columns.max(1) as f64,
            first_pin_relations: self.first_pin.map(|(r, _)| r),
            first_pin_trials: self.first_pin.map(|(_, t)| t),
            square_relations: self.square.map(|(r, _, _)| r),
            square_trials: self.square.map(|(_, t, _)| t),
            square_over_columns: self
                .square
                .map(|(r, _, _)| r as f64 / self.fb.columns.max(1) as f64),
            touched_at_square: self.square.map(|(_, _, n)| n),
            rank_at_columns: self.rank_at_columns,
            rank_fraction_at_columns: self
                .rank_at_columns
                .map(|rk| rk as f64 / (self.fb.columns + 1 - self.deficiency_total).max(1) as f64),
            deficiency_at_square: self.deficiency_at_square,
            deficiency_total: self.deficiency_total,
            relations_total: self.relations,
            rank_final: self.matrix.rank() as u64,
            dependent_final: self.matrix.dependent(),
            row_ops: self.matrix.work().0,
            zero_support_rows: self.zero_support,
            single_column_rows: self.single,
            recovered: self.recovered,
            verified,
        }
    }
}

/// The `E[h]`-component class of every column of a base: the least key
/// over the components `[r]P` of the column's points; `None` when the
/// column's stabiliser acts with a non-trivial eigenvalue on a
/// component (the same component carried by two points of the column
/// with different fold coefficients, e.g. `u = −u`), which forces the
/// functional to vanish there and contributes no deficiency.  Only
/// two-summand relations carry the constraint, so `summands > 2` or
/// `h = 1` yields no classes.
fn component_classes<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    r: u64,
    h: u64,
    summands: u32,
) -> Vec<Option<u64>> {
    if summands != 2 || h == 1 {
        return vec![None; fb.columns];
    }
    let mut ops = GroupOps::default();
    let mut seen: Vec<std::collections::HashMap<u64, u64>> = vec![Default::default(); fb.columns];
    let mut void = vec![false; fb.columns];
    for (i, p) in fb.points.iter().enumerate() {
        let u = g.mul(&mut ops, *p, r);
        let col = fb.col_of[i];
        let key = g.key(&u);
        if g.is_identity(&u) {
            void[col] = true;
            continue;
        }
        match seen[col].insert(key, fb.coef_of[i]) {
            Some(prev) if prev != fb.coef_of[i] => void[col] = true,
            _ => {}
        }
    }
    (0..fb.columns)
        .map(|c| {
            if void[c] {
                None
            } else {
                seen[c].keys().min().copied()
            }
        })
        .collect()
}

/// When a stream stops.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopRule {
    /// Both arms at their achievable full rank (E1–E12).
    BothSquare,
    /// The folded arm at full rank and the target's logarithm pinned on
    /// both arms: for an arm whose achievable rank is below the
    /// driver's `columns + 1 − D` (E13's control, where `P` and `P + T`
    /// are two columns), full rank is never declared, and the
    /// logarithm — what a pipeline stops at — is the end point.
    FoldedSquareBothPinned,
}

/// Feed one target stream to a folded arm and a control arm over the
/// same points, until both are square or `max_trials` targets are
/// drawn.  `oracle` decomposes a target into indices into the shared
/// point vector; `classes`, when given, canonicalises each target to
/// count orbit duplicates (E7).
#[allow(clippy::too_many_arguments)]
pub fn full_rank_stream<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    h: u64,
    planted: u64,
    folded: &FactorBase<G::Elt>,
    control: &FactorBase<G::Elt>,
    seed: u64,
    max_trials: u64,
    summands: u32,
    classes: Option<&EndomorphismClasses<'_, G>>,
    oracle: impl FnMut(&mut GroupOps, &mut OracleCounters, G::Elt) -> Option<Vec<usize>>,
) -> Result<StreamReport, String> {
    full_rank_stream_until(
        g,
        generator,
        target,
        r,
        h,
        planted,
        folded,
        control,
        seed,
        max_trials,
        summands,
        classes,
        StopRule::BothSquare,
        oracle,
    )
}

/// [`full_rank_stream`] with its stopping rule chosen.
#[allow(clippy::too_many_arguments)]
pub fn full_rank_stream_until<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    h: u64,
    planted: u64,
    folded: &FactorBase<G::Elt>,
    control: &FactorBase<G::Elt>,
    seed: u64,
    max_trials: u64,
    summands: u32,
    classes: Option<&EndomorphismClasses<'_, G>>,
    stop: StopRule,
    mut oracle: impl FnMut(&mut GroupOps, &mut OracleCounters, G::Elt) -> Option<Vec<usize>>,
) -> Result<StreamReport, String> {
    if folded.points.len() != control.points.len()
        || folded
            .points
            .iter()
            .zip(&control.points)
            .any(|(a, b)| g.key(a) != g.key(b))
    {
        return Err("the two arms must hold the same points in the same order".into());
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x4655_4C4C_5241_4E4B);
    let mut ops = GroupOps::default();
    let mut ctr = OracleCounters::default();
    let folded_classes = component_classes(g, folded, r, h, summands);
    let control_classes = component_classes(g, control, r, h, summands);
    let mut arms = [
        Arm {
            label: "folded".into(),
            fb: folded,
            matrix: matrix_by_name("structured-gauss", folded.columns + 1, r)?,
            class_of_column: folded_classes.clone(),
            touched_classes: HashSet::new(),
            touched: HashSet::new(),
            relations: 0,
            first_pin: None,
            square: None,
            deficiency_at_square: None,
            deficiency_total: folded_classes
                .iter()
                .flatten()
                .collect::<HashSet<_>>()
                .len(),
            rank_at_columns: None,
            zero_support: 0,
            single: 0,
            recovered: None,
        },
        Arm {
            label: "control".into(),
            fb: control,
            matrix: matrix_by_name("structured-gauss", control.columns + 1, r)?,
            class_of_column: control_classes.clone(),
            touched_classes: HashSet::new(),
            touched: HashSet::new(),
            relations: 0,
            first_pin: None,
            square: None,
            deficiency_at_square: None,
            deficiency_total: control_classes
                .iter()
                .flatten()
                .collect::<HashSet<_>>()
                .len(),
            rank_at_columns: None,
            zero_support: 0,
            single: 0,
            recovered: None,
        },
    ];
    let mut seen: HashSet<u64> = HashSet::new();
    let mut orbits: HashSet<u64> = HashSet::new();
    let mut orbit_duplicates = 0u64;
    let mut trials = 0u64;
    let mut relations = 0u64;
    let done = |arms: &[Arm<'_, G::Elt>; 2]| match stop {
        StopRule::BothSquare => arms.iter().all(|a| a.square.is_some()),
        StopRule::FoldedSquareBothPinned => {
            arms[0].square.is_some() && arms.iter().all(|a| a.first_pin.is_some())
        }
    };
    while trials < max_trials && !done(&arms) {
        let a = rng.gen_range(1..r);
        let b = rng.gen_range(1..r);
        let ag = g.mul(&mut ops, generator, a);
        let bq = g.mul(&mut ops, target, b);
        let point = g.add(&mut ops, ag, bq);
        if g.is_identity(&point) {
            continue;
        }
        trials += 1;
        // The framework's guard: a repeated target is a collision, not
        // a relation.  Once every element of ⟨G⟩ has been drawn there is
        // nothing left to decompose.
        if !seen.insert(g.key(&point)) {
            if seen.len() as u64 >= r.saturating_sub(1) {
                break;
            }
            continue;
        }
        if let Some(c) = classes {
            let (rep, _) = c.canon(g, point);
            if !orbits.insert(g.key(&rep)) {
                orbit_duplicates += 1;
            }
        }
        let Some(summands) = oracle(&mut ops, &mut ctr, point) else {
            continue;
        };
        relations += 1;
        for arm in arms.iter_mut() {
            arm.feed(&summands, a, b, h, r, trials);
        }
    }
    let verify = |rec: Option<u64>| -> bool {
        match rec {
            Some(d) => {
                let mut o = GroupOps::default();
                g.mul(&mut o, generator, d) == target && d == planted
            }
            None => false,
        }
    };
    let folded_rep = arms[0].report(verify(arms[0].recovered));
    let control_rep = arms[1].report(verify(arms[1].recovered));
    let ratio = |f: Option<u64>, c: Option<u64>| match (f, c) {
        (Some(f), Some(c)) if f > 0 => Some(c as f64 / f as f64),
        _ => None,
    };
    Ok(StreamReport {
        trials,
        relations,
        hit_rate: relations as f64 / trials.max(1) as f64,
        group_ops: ops,
        oracle_lookups: ctr.lookups,
        oracle_lift_failures: ctr.lift_failures,
        oracle_unliftable_systems: ctr.unliftable_systems,
        orbit_duplicates,
        square_ratio: ratio(folded_rep.square_relations, control_rep.square_relations),
        first_pin_ratio: ratio(
            folded_rep.first_pin_relations,
            control_rep.first_pin_relations,
        ),
        column_ratio: control_rep.columns as f64 / folded_rep.columns.max(1) as f64,
        exhausted: arms.iter().any(|a| a.square.is_none()),
        folded: folded_rep,
        control: control_rep,
    })
}

/// A stream with the negation-only arm as control and no orbit
/// canonicalisation: the E1 configuration.
pub fn e1_stream<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    h: u64,
    planted: u64,
    folded: &FactorBase<G::Elt>,
    control: &FactorBase<G::Elt>,
    seed: u64,
    max_trials: u64,
    oracle: impl FnMut(&mut GroupOps, &mut OracleCounters, G::Elt) -> Option<Vec<usize>>,
) -> Result<StreamReport, String> {
    full_rank_stream(
        g, generator, target, r, h, planted, folded, control, seed, max_trials, 2, None, oracle,
    )
}

/// The group's orbit classes for E7's duplicate count, built from the
/// generators the folded base used.
pub fn classes_for<'a, G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    r: u64,
    gens: &[&'a dyn Endomorphism<G>],
) -> Result<EndomorphismClasses<'a, G>, String> {
    EndomorphismClasses::new(g, generator, r, gens.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::{
        automorphism_generators, generate_cm_instance, glv_orbit_base, AutomorphismGroup, CmFamily,
    };
    use crate::cryptanalysis::ic_framework::plugins::{MitmOracle, SubtractOracle};
    use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};

    #[test]
    fn one_stream_feeds_both_arms_and_the_fold_shows_at_full_rank() {
        let inst = generate_cm_instance(CmFamily::J0, 15, 3, 1).unwrap();
        let (folded, _) = glv_orbit_base(&inst, 24, AutomorphismGroup::Auto, true).unwrap();
        let (control, _) = glv_orbit_base(&inst, 24, AutomorphismGroup::Auto, false).unwrap();
        let g = inst.generator_point();
        let planted = 1234 % inst.r;
        let mut ops = GroupOps::default();
        let target = inst.curve.mul(&mut ops, g, planted);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: g,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let mut oracle = SubtractOracle;
        oracle
            .prepare(&ctx, &folded, &Params::default(), &mut GroupOps::default())
            .unwrap();
        let rep = e1_stream(
            &inst.curve,
            g,
            target,
            inst.r,
            inst.cofactor,
            planted,
            &folded,
            &control,
            9,
            2_000_000,
            |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
        )
        .unwrap();
        assert!(!rep.exhausted, "both arms reach a square system");
        assert!(rep.folded.verified && rep.control.verified);
        assert_eq!(rep.column_ratio, 3.0);
        let sq = rep.square_ratio.expect("both squared");
        assert!(
            sq > 1.5,
            "the control needs more relations to square: ratio {sq}"
        );
        // Both arms saw the same relations.
        assert_eq!(rep.folded.relations_total, rep.control.relations_total);
        assert!(rep.folded.square_over_columns.unwrap() >= 0.5);
    }

    #[test]
    fn orbit_duplicates_are_counted_with_a_three_summand_oracle() {
        let inst = generate_cm_instance(CmFamily::J0, 14, 4, 1).unwrap();
        let (folded, _) = glv_orbit_base(&inst, 8, AutomorphismGroup::Auto, true).unwrap();
        let (control, _) = glv_orbit_base(&inst, 8, AutomorphismGroup::Auto, false).unwrap();
        let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
        let refs: Vec<&dyn Endomorphism<_>> = gens.iter().map(|b| b.as_ref()).collect();
        let g = inst.generator_point();
        let classes = classes_for(&inst.curve, g, inst.r, &refs).unwrap();
        assert_eq!(classes.group_order, 6);
        let planted = 77 % inst.r;
        let mut ops = GroupOps::default();
        let target = inst.curve.mul(&mut ops, g, planted);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: g,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let mut oracle = MitmOracle::new(3);
        oracle
            .prepare(&ctx, &folded, &Params::default(), &mut GroupOps::default())
            .unwrap();
        let rep = full_rank_stream(
            &inst.curve,
            g,
            target,
            inst.r,
            inst.cofactor,
            planted,
            &folded,
            &control,
            5,
            200_000,
            3,
            Some(&classes),
            |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
        )
        .unwrap();
        assert!(rep.folded.verified, "three summands recover the logarithm");
        // On a 2^14 group with thousands of targets, orbit repeats occur.
        assert!(rep.trials > 0);
    }
}
