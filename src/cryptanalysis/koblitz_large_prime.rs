//! Exact, full-width elimination of one- and two-large-prime relations.
//!
//! A candidate first proves its point equation in the prime-order subgroup.
//! Each residual point is then represented by its signed Frobenius orbit and
//! an exact coefficient modulo the subgroup order. Pivot reduction preserves
//! the coefficients of the original public target and generator. A completed
//! row is checked again as a group equation before it can enter linear
//! algebra. This module does not discover partial relations or time a DLP.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::binary_ecc::curve::point_neg;
use crate::binary_ecc::{BinaryPoint, F2mElement};
use crate::cryptanalysis::koblitz_index_calculus::{point_key, FrobeniusFactorBase, KoblitzCurve};
use crate::utils::mod_inverse;

type PointKey = (BigUint, BigUint);

/// A public, independently checkable partial point equation:
/// [a]G + [b]Q = sum(small points) + sum(residual points).
#[derive(Clone, Debug)]
pub struct PartialCandidate {
    pub source_id: u64,
    pub coef_a: BigUint,
    pub coef_b: BigUint,
    pub small_point_indices: Vec<usize>,
    pub residual_points: Vec<BinaryPoint>,
}

/// Bounds on the graph rather than silent truncations of its algebra.
#[derive(Clone, Copy, Debug)]
pub struct LargePrimeLimits {
    pub max_merge_depth: usize,
    pub max_row_terms: usize,
    pub max_pivots: usize,
}

impl Default for LargePrimeLimits {
    fn default() -> Self {
        Self {
            max_merge_depth: 32,
            max_row_terms: 4096,
            max_pivots: 1_000_000,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CapReason {
    MergeDepth,
    RowTerms,
    PivotCount,
}

/// The graph makes no completion claim when a configured cap is reached.
#[derive(Clone, Debug)]
pub enum FeedOutcome {
    PivotStored {
        merge_depth: usize,
    },
    Full(WideFullRelation),
    UnknownCap {
        reason: CapReason,
        merge_depth: usize,
    },
}

/// A row over signed factor-base orbit representatives only.
///
/// Source coefficients make the row independently replayable from retained
/// input partials. The small coefficients and both right-hand-side scalars
/// are always reduced modulo the full subgroup order.
#[derive(Clone, Debug)]
pub struct WideFullRelation {
    pub coef_a: BigUint,
    pub coef_b: BigUint,
    pub small: BTreeMap<usize, BigUint>,
    pub source_coeffs: BTreeMap<u64, BigUint>,
    pub merge_depth: usize,
}

impl WideFullRelation {
    pub fn dense_row(&self, columns: usize) -> Result<Vec<BigUint>, &'static str> {
        let mut row = vec![BigUint::zero(); columns];
        for (&column, coefficient) in &self.small {
            let entry = row.get_mut(column).ok_or("small column out of range")?;
            *entry = coefficient.clone();
        }
        Ok(row)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct LargePrimeStats {
    pub accepted_candidates: u64,
    pub pivot_rows: usize,
    pub full_rows: u64,
    pub full_from_partials: u64,
    pub capped_rows: u64,
    pub pivot_reductions: u64,
}

#[derive(Clone, Debug)]
struct SparseRow {
    coef_a: BigUint,
    coef_b: BigUint,
    small: BTreeMap<usize, BigUint>,
    large: BTreeMap<PointKey, BigUint>,
    sources: BTreeMap<u64, BigUint>,
}

fn add_mod<K: Ord>(
    map: &mut BTreeMap<K, BigUint>,
    key: K,
    coefficient: &BigUint,
    modulus: &BigUint,
) {
    let value = (map.get(&key).cloned().unwrap_or_default() + coefficient) % modulus;
    if value.is_zero() {
        map.remove(&key);
    } else {
        map.insert(key, value);
    }
}

fn sub_mod(a: &BigUint, b: &BigUint, modulus: &BigUint) -> BigUint {
    (a + modulus - b) % modulus
}

impl SparseRow {
    fn term_count(&self) -> usize {
        self.small.len() + self.large.len() + self.sources.len()
    }

    fn scale(&mut self, factor: &BigUint, modulus: &BigUint) {
        self.coef_a = (&self.coef_a * factor) % modulus;
        self.coef_b = (&self.coef_b * factor) % modulus;
        for value in self.small.values_mut() {
            *value = (&*value * factor) % modulus;
        }
        for value in self.sources.values_mut() {
            *value = (&*value * factor) % modulus;
        }
        for value in self.large.values_mut() {
            *value = (&*value * factor) % modulus;
        }
    }

    fn subtract_scaled(&mut self, pivot: &SparseRow, factor: &BigUint, modulus: &BigUint) {
        self.coef_a = sub_mod(&self.coef_a, &((&pivot.coef_a * factor) % modulus), modulus);
        self.coef_b = sub_mod(&self.coef_b, &((&pivot.coef_b * factor) % modulus), modulus);
        for (&key, value) in &pivot.small {
            let neg = sub_mod(&BigUint::zero(), &((value * factor) % modulus), modulus);
            add_mod(&mut self.small, key, &neg, modulus);
        }
        for (key, value) in &pivot.large {
            let neg = sub_mod(&BigUint::zero(), &((value * factor) % modulus), modulus);
            add_mod(&mut self.large, key.clone(), &neg, modulus);
        }
        for (&key, value) in &pivot.sources {
            let neg = sub_mod(&BigUint::zero(), &((value * factor) % modulus), modulus);
            add_mod(&mut self.sources, key, &neg, modulus);
        }
    }
}

enum Reduced {
    Pivot(PointKey, usize),
    Full(SparseRow, usize),
    Capped(CapReason, usize),
}

struct Reducer {
    modulus: BigUint,
    limits: LargePrimeLimits,
    pivots: BTreeMap<PointKey, SparseRow>,
    reductions: u64,
}

impl Reducer {
    fn new(modulus: BigUint, limits: LargePrimeLimits) -> Result<Self, &'static str> {
        if modulus < BigUint::from(2u32)
            || limits.max_merge_depth == 0
            || limits.max_row_terms == 0
            || limits.max_pivots == 0
        {
            return Err("invalid large-prime modulus or cap");
        }
        Ok(Self {
            modulus,
            limits,
            pivots: BTreeMap::new(),
            reductions: 0,
        })
    }

    fn feed(&mut self, mut row: SparseRow) -> Result<Reduced, &'static str> {
        if row.term_count() > self.limits.max_row_terms {
            return Ok(Reduced::Capped(CapReason::RowTerms, 0));
        }
        let mut depth = 0;
        loop {
            if row.large.is_empty() {
                return Ok(Reduced::Full(row, depth));
            }
            let existing = row
                .large
                .iter()
                .find(|(key, _)| self.pivots.contains_key(*key))
                .map(|(key, value)| (key.clone(), value.clone()));
            if let Some((key, factor)) = existing {
                if depth >= self.limits.max_merge_depth {
                    return Ok(Reduced::Capped(CapReason::MergeDepth, depth));
                }
                row.subtract_scaled(&self.pivots[&key], &factor, &self.modulus);
                depth += 1;
                self.reductions += 1;
                if row.term_count() > self.limits.max_row_terms {
                    return Ok(Reduced::Capped(CapReason::RowTerms, depth));
                }
                continue;
            }
            if self.pivots.len() >= self.limits.max_pivots {
                return Ok(Reduced::Capped(CapReason::PivotCount, depth));
            }
            let (key, coefficient) = row.large.iter().next().expect("nonempty row");
            let key = key.clone();
            let inverse = mod_inverse(coefficient, &self.modulus)
                .ok_or("large-prime pivot is not invertible modulo subgroup order")?;
            row.scale(&inverse, &self.modulus);
            self.pivots.insert(key.clone(), row);
            return Ok(Reduced::Pivot(key, depth));
        }
    }
}

/// A wide, signed-Frobenius large-prime adapter for subgroup point relations.
///
/// Its input is a completed decomposition with zero, one or two residual
/// points. A discovery solver and a source-pinned total-runtime harness are
/// separate gates.
pub struct WideLargePrimeEliminator<'a> {
    curve: &'a KoblitzCurve,
    base: &'a FrobeniusFactorBase,
    target: BinaryPoint,
    base_index: HashMap<PointKey, usize>,
    reducer: Reducer,
    seen_sources: BTreeSet<u64>,
    stats: LargePrimeStats,
}

impl<'a> WideLargePrimeEliminator<'a> {
    pub fn new(
        curve: &'a KoblitzCurve,
        base: &'a FrobeniusFactorBase,
        target: BinaryPoint,
        limits: LargePrimeLimits,
    ) -> Result<Self, &'static str> {
        let modulus = &curve.subgroup_order;
        if !curve.frobenius_is_endomorphism
            || curve.lambda.modpow(&BigUint::from(curve.n), modulus) != BigUint::one()
            || target == BinaryPoint::Infinity
            || !curve.curve.is_on_curve(&target)
            || curve.mul(&target, modulus) != BinaryPoint::Infinity
            || !curve.curve.is_on_curve(curve.generator())
            || curve.mul(curve.generator(), modulus) != BinaryPoint::Infinity
            || base.signed_orbit_of.len() != base.points.len()
            || base.signed_orbits.is_empty()
        {
            return Err("large-prime adapter requires a valid subgroup and orbit base");
        }
        for orbit in &base.signed_orbits {
            let index = *orbit.first().ok_or("empty signed orbit")?;
            let rep = base.points.get(index).ok_or("invalid signed orbit index")?;
            if *rep == BinaryPoint::Infinity
                || !curve.curve.is_on_curve(rep)
                || curve.mul(rep, modulus) != BinaryPoint::Infinity
            {
                return Err("signed orbit representative is outside the subgroup");
            }
        }
        Ok(Self {
            curve,
            base,
            target,
            base_index: base.index_map(),
            reducer: Reducer::new(modulus.clone(), limits)?,
            seen_sources: BTreeSet::new(),
            stats: LargePrimeStats::default(),
        })
    }

    pub fn stats(&self) -> LargePrimeStats {
        let mut stats = self.stats;
        stats.pivot_rows = self.reducer.pivots.len();
        stats.pivot_reductions = self.reducer.reductions;
        stats
    }

    fn canonical_residual(&self, point: &BinaryPoint) -> Result<(PointKey, BigUint), &'static str> {
        let modulus = &self.reducer.modulus;
        let mut walk = point.clone();
        let mut best: Option<(PointKey, BinaryPoint, u32, bool)> = None;
        for phase in 0..self.curve.n {
            for negated in [false, true] {
                let candidate = if negated {
                    point_neg(&walk)
                } else {
                    walk.clone()
                };
                let key = point_key(&candidate);
                if best.as_ref().is_none_or(|(old, _, _, _)| key < *old) {
                    best = Some((key, candidate, phase, negated));
                }
            }
            walk = self.curve.frobenius(&walk);
        }
        let (key, rep, phase, negated) = best.ok_or("empty Frobenius orbit")?;
        let phase_factor = self.curve.lambda.modpow(&BigUint::from(phase), modulus);
        let inverse =
            mod_inverse(&phase_factor, modulus).ok_or("Frobenius phase is not invertible")?;
        let coefficient = if negated {
            sub_mod(&BigUint::zero(), &inverse, modulus)
        } else {
            inverse
        };
        if self.curve.mul(&rep, &coefficient) != *point {
            return Err("residual Frobenius lifting failed");
        }
        Ok((key, coefficient))
    }

    fn point_from_key(&self, key: &PointKey) -> Result<BinaryPoint, &'static str> {
        if key.0.is_zero() {
            return Err("identity is not a large-prime column");
        }
        Ok(BinaryPoint::Affine {
            x: F2mElement::from_biguint(&(&key.0 - BigUint::one()), self.curve.n),
            y: F2mElement::from_biguint(&key.1, self.curve.n),
        })
    }

    fn verify_row(&self, row: &SparseRow) -> Result<(), &'static str> {
        let mut lhs = self.curve.mul(self.curve.generator(), &row.coef_a);
        lhs = self
            .curve
            .add(&lhs, &self.curve.mul(&self.target, &row.coef_b));
        let mut rhs = BinaryPoint::Infinity;
        for (&orbit, coefficient) in &row.small {
            let index = *self
                .base
                .signed_orbits
                .get(orbit)
                .and_then(|points| points.first())
                .ok_or("small orbit index out of range")?;
            rhs = self
                .curve
                .add(&rhs, &self.curve.mul(&self.base.points[index], coefficient));
        }
        for (key, coefficient) in &row.large {
            let rep = self.point_from_key(key)?;
            rhs = self.curve.add(&rhs, &self.curve.mul(&rep, coefficient));
        }
        if lhs != rhs {
            return Err("large-prime row failed group replay");
        }
        Ok(())
    }

    pub fn feed(&mut self, candidate: PartialCandidate) -> Result<FeedOutcome, &'static str> {
        if self.seen_sources.contains(&candidate.source_id) {
            return Err("duplicate partial source id");
        }
        if candidate.residual_points.len() > 2 {
            return Err("at most two residual points are supported");
        }
        let modulus = &self.reducer.modulus;
        let coef_a = &candidate.coef_a % modulus;
        let coef_b = &candidate.coef_b % modulus;
        let mut rhs = BinaryPoint::Infinity;
        let mut small = BTreeMap::new();
        for &index in &candidate.small_point_indices {
            let point = self
                .base
                .points
                .get(index)
                .ok_or("small point index out of range")?;
            if *point == BinaryPoint::Infinity
                || !self.curve.curve.is_on_curve(point)
                || self.curve.mul(point, modulus) != BinaryPoint::Infinity
            {
                return Err("small point is outside the subgroup");
            }
            rhs = self.curve.add(&rhs, point);
            let (orbit, phase, negated) = *self
                .base
                .signed_orbit_of
                .get(index)
                .ok_or("small point lacks signed orbit label")?;
            let rep_index = *self
                .base
                .signed_orbits
                .get(orbit)
                .and_then(|points| points.first())
                .ok_or("small point has invalid signed orbit label")?;
            let phase_factor = self.curve.lambda.modpow(&BigUint::from(phase), modulus);
            let coefficient = if negated {
                sub_mod(&BigUint::zero(), &phase_factor, modulus)
            } else {
                phase_factor
            };
            if self.curve.mul(&self.base.points[rep_index], &coefficient) != *point {
                return Err("small signed orbit label failed group replay");
            }
            add_mod(&mut small, orbit, &coefficient, modulus);
        }
        let mut large = BTreeMap::new();
        for point in &candidate.residual_points {
            if *point == BinaryPoint::Infinity
                || !self.curve.curve.is_on_curve(point)
                || self.curve.mul(point, modulus) != BinaryPoint::Infinity
                || self.base_index.contains_key(&point_key(point))
            {
                return Err("residual point is invalid or already in the factor base");
            }
            rhs = self.curve.add(&rhs, point);
            let (key, coefficient) = self.canonical_residual(point)?;
            if self.base_index.contains_key(&key) {
                return Err("residual orbit intersects the factor base");
            }
            add_mod(&mut large, key, &coefficient, modulus);
        }
        let lhs = self.curve.add(
            &self.curve.mul(self.curve.generator(), &coef_a),
            &self.curve.mul(&self.target, &coef_b),
        );
        if lhs != rhs {
            return Err("partial relation failed exact point-sum replay");
        }
        let mut sources = BTreeMap::new();
        sources.insert(candidate.source_id, BigUint::one());
        let row = SparseRow {
            coef_a,
            coef_b,
            small,
            large,
            sources,
        };
        self.verify_row(&row)?;
        let reduced = self.reducer.feed(row)?;
        let outcome = match reduced {
            Reduced::Pivot(key, merge_depth) => {
                if let Err(error) = self.verify_row(&self.reducer.pivots[&key]) {
                    self.reducer.pivots.remove(&key);
                    return Err(error);
                }
                FeedOutcome::PivotStored { merge_depth }
            }
            Reduced::Full(row, merge_depth) => {
                self.verify_row(&row)?;
                self.stats.full_rows += 1;
                if merge_depth > 0 {
                    self.stats.full_from_partials += 1;
                }
                FeedOutcome::Full(WideFullRelation {
                    coef_a: row.coef_a,
                    coef_b: row.coef_b,
                    small: row.small,
                    source_coeffs: row.sources,
                    merge_depth,
                })
            }
            Reduced::Capped(reason, merge_depth) => {
                self.stats.capped_rows += 1;
                FeedOutcome::UnknownCap {
                    reason,
                    merge_depth,
                }
            }
        };
        self.seen_sources.insert(candidate.source_id);
        self.stats.accepted_candidates += 1;
        Ok(outcome)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_index_calculus::build_explicit_frobenius_orbit_factor_base;

    fn fake_key(x: u32) -> PointKey {
        (BigUint::from(x), BigUint::zero())
    }

    fn algebra_row(id: u64, a: u64, small: usize, large: &[u32]) -> SparseRow {
        SparseRow {
            coef_a: BigUint::from(a),
            coef_b: BigUint::zero(),
            small: BTreeMap::from([(small, BigUint::one())]),
            large: large
                .iter()
                .map(|&x| (fake_key(x), BigUint::one()))
                .collect(),
            sources: BTreeMap::from([(id, BigUint::one())]),
        }
    }

    #[test]
    fn full_width_double_large_prime_cycle_has_exact_source_coefficients() {
        let modulus: BigUint = (BigUint::one() << 81usize) - BigUint::one();
        let mut reducer = Reducer::new(modulus.clone(), LargePrimeLimits::default()).unwrap();
        assert!(matches!(
            reducer.feed(algebra_row(1, 10, 0, &[2, 3])).unwrap(),
            Reduced::Pivot(_, 0)
        ));
        assert!(matches!(
            reducer.feed(algebra_row(2, 4, 1, &[2])).unwrap(),
            Reduced::Pivot(_, 1)
        ));
        let Reduced::Full(row, depth) = reducer.feed(algebra_row(3, 7, 2, &[3])).unwrap() else {
            panic!("three partials must close the two-residual cycle");
        };
        assert!(depth > 0);
        assert!(row.large.is_empty());
        assert_eq!(row.sources.len(), 3);
        let a_from_sources = row
            .sources
            .iter()
            .fold(BigUint::zero(), |acc, (id, coeff)| {
                let input_a = match id {
                    1 => 10u32,
                    2 => 4,
                    3 => 7,
                    _ => panic!("unknown source"),
                };
                (acc + coeff * BigUint::from(input_a)) % &modulus
            });
        assert_eq!(a_from_sources, row.coef_a);
        assert!(row.coef_a.bits() <= 81);
    }

    #[test]
    fn graph_caps_are_unknown_and_do_not_store_unreduced_rows() {
        let modulus = BigUint::from(101u32);
        let limits = LargePrimeLimits {
            max_merge_depth: 1,
            max_row_terms: 20,
            max_pivots: 2,
        };
        let mut reducer = Reducer::new(modulus.clone(), limits).unwrap();
        assert!(matches!(
            reducer.feed(algebra_row(1, 1, 0, &[2])).unwrap(),
            Reduced::Pivot(_, 0)
        ));
        assert!(matches!(
            reducer.feed(algebra_row(2, 2, 1, &[3])).unwrap(),
            Reduced::Pivot(_, 0)
        ));
        assert!(matches!(
            reducer.feed(algebra_row(3, 3, 2, &[2, 3])).unwrap(),
            Reduced::Capped(CapReason::MergeDepth, 1)
        ));
        assert_eq!(reducer.pivots.len(), 2);
        assert!(matches!(
            reducer.feed(algebra_row(4, 4, 3, &[4])).unwrap(),
            Reduced::Capped(CapReason::PivotCount, 0)
        ));
        assert_eq!(reducer.pivots.len(), 2);
        let mut narrow = Reducer::new(
            modulus,
            LargePrimeLimits {
                max_row_terms: 2,
                ..LargePrimeLimits::default()
            },
        )
        .unwrap();
        assert!(matches!(
            narrow.feed(algebra_row(5, 5, 4, &[5])).unwrap(),
            Reduced::Capped(CapReason::RowTerms, 0)
        ));
        assert!(narrow.pivots.is_empty());
    }

    #[test]
    fn group_replay_closes_two_residual_orbits_and_rejects_bad_inputs() {
        let curve = KoblitzCurve::new(1, 11).unwrap();
        let BinaryPoint::Affine { x, .. } = curve.generator() else {
            panic!("generator is affine");
        };
        let base = build_explicit_frobenius_orbit_factor_base(&curve, &[x.clone()]).unwrap();
        let target = curve.mul(curve.generator(), &BigUint::from(17u32));
        let mut graph =
            WideLargePrimeEliminator::new(&curve, &base, target, LargePrimeLimits::default())
                .unwrap();
        let mut residuals: Vec<(u32, BinaryPoint, PointKey)> = Vec::new();
        for scalar in 2..200u32 {
            let point = curve.mul(curve.generator(), &BigUint::from(scalar));
            if base.index_map().contains_key(&point_key(&point)) {
                continue;
            }
            let (key, _) = graph.canonical_residual(&point).unwrap();
            if residuals.iter().all(|(_, _, existing)| *existing != key) {
                residuals.push((scalar, point, key));
            }
            if residuals.len() == 2 {
                break;
            }
        }
        assert_eq!(residuals.len(), 2);
        let r = &curve.subgroup_order;
        let (key, positive) = graph.canonical_residual(&residuals[0].1).unwrap();
        let (negative_key, negative) = graph
            .canonical_residual(&point_neg(&residuals[0].1))
            .unwrap();
        assert_eq!(key, negative_key);
        assert_eq!((positive + negative) % r, BigUint::zero());
        let make = |source_id: u64, b: u32, selected: &[usize]| {
            let sum = selected.iter().fold(BigUint::one(), |acc, &i| {
                acc + BigUint::from(residuals[i].0)
            });
            PartialCandidate {
                source_id,
                coef_a: (sum + r - (BigUint::from(b) * BigUint::from(17u32)) % r) % r,
                coef_b: BigUint::from(b),
                small_point_indices: vec![base.index_map()[&point_key(curve.generator())]],
                residual_points: selected.iter().map(|&i| residuals[i].1.clone()).collect(),
            }
        };
        let mut bad = make(99, 1, &[0]);
        bad.coef_a += BigUint::one();
        assert!(graph.feed(bad).is_err());
        assert!(matches!(
            graph.feed(make(1, 1, &[0, 1])).unwrap(),
            FeedOutcome::PivotStored { .. }
        ));
        assert!(matches!(
            graph.feed(make(2, 2, &[0])).unwrap(),
            FeedOutcome::PivotStored { .. }
        ));
        let FeedOutcome::Full(full) = graph.feed(make(3, 3, &[1])).unwrap() else {
            panic!("two-residual graph did not close");
        };
        assert_eq!(full.source_coeffs.len(), 3);
        assert!(full.merge_depth > 0);
        assert_eq!(
            full.dense_row(base.unknowns()).unwrap().len(),
            base.unknowns()
        );
        let stats = graph.stats();
        assert_eq!(stats.accepted_candidates, 3);
        assert_eq!(stats.full_from_partials, 1);
        assert_eq!(stats.full_rows, 1);
        assert_eq!(stats.pivot_rows, 2);
        assert!(graph.feed(make(3, 3, &[1])).is_err(), "duplicate source id");
    }
}
