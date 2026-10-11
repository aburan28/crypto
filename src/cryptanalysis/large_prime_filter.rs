//! Reusable large-prime graph and hypergraph filtering for ECDLP relations.
//!
//! Collectors represent a partial relation as a sparse equation whose
//! low-numbered columns are factor-base logarithms and whose `large` columns
//! are temporary large primes.  This module removes the temporary columns in
//! two exact stages:
//!
//! 1. peel the degree-one vertices of the large-prime hypergraph; and
//! 2. run sparse, Markowitz-ordered elimination on its two-core.
//!
//! Every emitted complete relation carries a sparse linear combination of the
//! input relation ids.  A caller can therefore replay the original curve-group
//! witnesses instead of trusting the filter.  The arithmetic supports one,
//! two, or arbitrarily many large primes per relation and arbitrary nonzero
//! coefficients over an odd prime field.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

#[inline]
fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 + b as u128) % modulus as u128) as u64
}

#[inline]
fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

#[inline]
fn neg_mod(a: u64, modulus: u64) -> u64 {
    if a == 0 {
        0
    } else {
        modulus - a
    }
}

fn inverse_mod(a: u64, modulus: u64) -> Option<u64> {
    let (mut old_r, mut r) = (a, modulus);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q as i128 * s);
    }
    (old_r == 1).then(|| old_s.rem_euclid(modulus as i128) as u64)
}

fn normalise<K: Ord + Copy>(terms: &mut Vec<(K, u64)>, modulus: u64) {
    let mut merged = BTreeMap::<K, u64>::new();
    for &(key, coefficient) in terms.iter() {
        let value = coefficient % modulus;
        if value != 0 {
            let slot = merged.entry(key).or_default();
            *slot = add_mod(*slot, value, modulus);
        }
    }
    terms.clear();
    terms.extend(merged.into_iter().filter(|(_, value)| *value != 0));
}

/// One relation before its large-prime columns have been eliminated.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PartialRelation {
    /// Stable collector id used in replay certificates.
    pub id: u64,
    /// `(factor-base column, coefficient)` pairs.
    pub small: Vec<(u32, u64)>,
    /// `(large-prime identity, coefficient)` pairs.
    pub large: Vec<(u64, u64)>,
    pub rhs: u64,
}

impl PartialRelation {
    fn normalise(&mut self, modulus: u64) {
        normalise(&mut self.small, modulus);
        normalise(&mut self.large, modulus);
        self.rhs %= modulus;
    }
}

/// A relation containing factor-base columns only.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompleteRelation {
    pub small: Vec<(u32, u64)>,
    pub rhs: u64,
    /// Coefficients of input relation ids whose linear combination produced
    /// this row.  This is the replay witness for curve-group verification.
    pub sources: Vec<(u64, u64)>,
}

/// Safety limits for hypergraph elimination.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct LargePrimeFilterOptions {
    pub modulus: u64,
    /// Refuse a merge whose factor-base row would exceed this weight.
    pub max_small_weight: usize,
    /// Refuse a merge whose replay certificate would exceed this weight.
    pub max_sources: usize,
}

impl Default for LargePrimeFilterOptions {
    fn default() -> Self {
        Self {
            modulus: 0,
            max_small_weight: 4096,
            max_sources: 4096,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct LargePrimeFilterReport {
    pub input_relations: usize,
    pub input_large_vertices: usize,
    pub duplicate_relations: usize,
    pub peeled_relations: usize,
    pub core_relations: usize,
    pub core_large_vertices: usize,
    pub pivots: usize,
    pub merges: usize,
    pub rejected_fill: usize,
    pub complete_relations: usize,
    /// Complete rows produced from a nontrivial graph/hypergraph cycle.
    pub cycles_closed: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LargePrimeFilterError {
    InvalidModulus,
    DuplicateId(u64),
    NonInvertiblePivot { large_prime: u64, coefficient: u64 },
}

#[derive(Clone, Debug)]
struct WorkRow {
    small: Vec<(u32, u64)>,
    large: Vec<(u64, u64)>,
    rhs: u64,
    sources: Vec<(u64, u64)>,
}

impl WorkRow {
    fn from_partial(row: PartialRelation, modulus: u64) -> Self {
        Self {
            small: row.small,
            large: row.large,
            rhs: row.rhs,
            sources: vec![(row.id, 1 % modulus)],
        }
    }

    fn coefficient(&self, large_prime: u64) -> u64 {
        self.large
            .binary_search_by_key(&large_prime, |&(key, _)| key)
            .ok()
            .map_or(0, |index| self.large[index].1)
    }

    /// `self + scale * other`.
    fn axpy(&self, scale: u64, other: &Self, modulus: u64) -> Self {
        fn combine<K: Ord + Copy>(
            left: &[(K, u64)],
            scale: u64,
            right: &[(K, u64)],
            modulus: u64,
        ) -> Vec<(K, u64)> {
            let mut out = BTreeMap::<K, u64>::new();
            for &(key, value) in left {
                out.insert(key, value);
            }
            for &(key, value) in right {
                let slot = out.entry(key).or_default();
                *slot = add_mod(*slot, mul_mod(scale, value, modulus), modulus);
            }
            out.into_iter().filter(|(_, value)| *value != 0).collect()
        }
        Self {
            small: combine(&self.small, scale, &other.small, modulus),
            large: combine(&self.large, scale, &other.large, modulus),
            rhs: add_mod(self.rhs, mul_mod(scale, other.rhs, modulus), modulus),
            sources: combine(&self.sources, scale, &other.sources, modulus),
        }
    }

    fn complete(self) -> CompleteRelation {
        debug_assert!(self.large.is_empty());
        CompleteRelation {
            small: self.small,
            rhs: self.rhs,
            sources: self.sources,
        }
    }
}

fn duplicate_key(row: &PartialRelation) -> (Vec<(u32, u64)>, Vec<(u64, u64)>, u64) {
    (row.small.clone(), row.large.clone(), row.rhs)
}

/// Peel and eliminate a batch of partial relations.
pub fn filter_large_primes(
    mut relations: Vec<PartialRelation>,
    options: &LargePrimeFilterOptions,
) -> Result<(Vec<CompleteRelation>, LargePrimeFilterReport), LargePrimeFilterError> {
    let modulus = options.modulus;
    // The extended-Euclid implementation below keeps signed coefficients in
    // `i128`; the ECDLP sparse stack uses the same explicit `< 2^63` bound.
    if !(3..(1u64 << 63)).contains(&modulus) || modulus.is_multiple_of(2) {
        return Err(LargePrimeFilterError::InvalidModulus);
    }
    let mut report = LargePrimeFilterReport {
        input_relations: relations.len(),
        ..LargePrimeFilterReport::default()
    };
    let mut ids = BTreeSet::new();
    for row in &mut relations {
        if !ids.insert(row.id) {
            return Err(LargePrimeFilterError::DuplicateId(row.id));
        }
        row.normalise(modulus);
    }
    let mut seen = BTreeSet::new();
    relations.retain(|row| {
        let fresh = seen.insert(duplicate_key(row));
        report.duplicate_relations += usize::from(!fresh);
        fresh
    });

    let mut incident = HashMap::<u64, Vec<usize>>::new();
    for (index, row) in relations.iter().enumerate() {
        for &(large_prime, _) in &row.large {
            incident.entry(large_prime).or_default().push(index);
        }
    }
    report.input_large_vertices = incident.len();

    // Hypergraph two-core peeling. A relation touched by a degree-one vertex
    // cannot participate in a dependency that cancels every large prime.
    let mut active = vec![true; relations.len()];
    let mut degree: HashMap<u64, usize> = incident
        .iter()
        .map(|(&large_prime, rows)| (large_prime, rows.len()))
        .collect();
    let mut queue: VecDeque<u64> = degree
        .iter()
        .filter_map(|(&large_prime, &d)| (d == 1).then_some(large_prime))
        .collect();
    while let Some(large_prime) = queue.pop_front() {
        if degree.get(&large_prime).copied().unwrap_or(0) != 1 {
            continue;
        }
        let Some(row_index) = incident[&large_prime]
            .iter()
            .copied()
            .find(|&index| active[index])
        else {
            continue;
        };
        active[row_index] = false;
        report.peeled_relations += 1;
        for &(other, _) in &relations[row_index].large {
            if let Some(d) = degree.get_mut(&other) {
                *d = d.saturating_sub(1);
                if *d == 1 {
                    queue.push_back(other);
                }
            }
        }
    }

    let mut complete = Vec::new();
    let mut core = Vec::new();
    for (index, row) in relations.into_iter().enumerate() {
        if row.large.is_empty() {
            complete.push(WorkRow::from_partial(row, modulus).complete());
        } else if active[index] {
            core.push(WorkRow::from_partial(row, modulus));
        }
    }
    report.core_relations = core.len();
    let mut core_degree = HashMap::<u64, usize>::new();
    for row in &core {
        for &(large_prime, _) in &row.large {
            *core_degree.entry(large_prime).or_default() += 1;
        }
    }
    report.core_large_vertices = core_degree.len();

    // Light rows first.  When several columns are available, choose the
    // Markowitz minimum `(degree - 1) * (row_weight - 1)`.
    core.sort_by_key(|row| (row.large.len(), row.small.len(), row.sources[0].0));
    let mut pivots = HashMap::<u64, WorkRow>::new();
    for mut row in core {
        let mut rejected_for_fill = false;
        loop {
            let pivotable = row
                .large
                .iter()
                .filter(|(large_prime, _)| pivots.contains_key(large_prime))
                .min_by_key(|(large_prime, _)| {
                    let degree = core_degree.get(large_prime).copied().unwrap_or(1);
                    degree.saturating_sub(1) * row.large.len().saturating_sub(1)
                })
                .map(|&(large_prime, coefficient)| (large_prime, coefficient));
            let Some((large_prime, coefficient)) = pivotable else {
                break;
            };
            let pivot = &pivots[&large_prime];
            let pivot_coefficient = pivot.coefficient(large_prime);
            let inverse = inverse_mod(pivot_coefficient, modulus).ok_or(
                LargePrimeFilterError::NonInvertiblePivot {
                    large_prime,
                    coefficient: pivot_coefficient,
                },
            )?;
            let scale = neg_mod(mul_mod(coefficient, inverse, modulus), modulus);
            let merged = row.axpy(scale, pivot, modulus);
            if merged.small.len() > options.max_small_weight
                || merged.sources.len() > options.max_sources
            {
                report.rejected_fill += 1;
                rejected_for_fill = true;
                break;
            }
            row = merged;
            report.merges += 1;
        }
        // The row still contains a column already owned by `pivots`.  It may
        // not replace that pivot: doing so would lose the earlier row and make
        // subsequent eliminations order-dependent.  Dropping it is the
        // conservative meaning of the caller's fill bound.
        if rejected_for_fill {
            continue;
        }
        if row.large.is_empty() {
            if row.sources.len() > 1 {
                report.cycles_closed += 1;
            }
            complete.push(row.complete());
            continue;
        }
        let (&(large_prime, _), _) = row
            .large
            .iter()
            .map(|term| {
                let degree = core_degree.get(&term.0).copied().unwrap_or(1);
                let score = degree.saturating_sub(1) * row.large.len().saturating_sub(1);
                (term, score)
            })
            .min_by_key(|(_, score)| *score)
            .expect("a partial row has a large-prime column");
        pivots.insert(large_prime, row);
        report.pivots += 1;
    }
    report.complete_relations = complete.len();
    Ok((complete, report))
}

/// Recompute a complete relation from its source certificate.
pub fn replay_complete_relation(
    complete: &CompleteRelation,
    inputs: &[PartialRelation],
    modulus: u64,
) -> Option<CompleteRelation> {
    let by_id: HashMap<u64, &PartialRelation> = inputs.iter().map(|row| (row.id, row)).collect();
    let mut accumulated = WorkRow {
        small: Vec::new(),
        large: Vec::new(),
        rhs: 0,
        sources: Vec::new(),
    };
    for &(id, coefficient) in &complete.sources {
        let mut source = (*by_id.get(&id)?).clone();
        source.normalise(modulus);
        let row = WorkRow::from_partial(source, modulus);
        accumulated = accumulated.axpy(coefficient, &row, modulus);
    }
    accumulated.large.is_empty().then(|| accumulated.complete())
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 101;

    fn row(id: u64, small: &[(u32, u64)], large: &[(u64, u64)], rhs: u64) -> PartialRelation {
        PartialRelation {
            id,
            small: small.to_vec(),
            large: large.to_vec(),
            rhs,
        }
    }

    #[test]
    fn single_large_prime_pair_closes_and_replays() {
        let inputs = vec![
            row(10, &[(0, 1)], &[(700, 1)], 8),
            row(11, &[(1, 1)], &[(700, 1)], 13),
        ];
        let options = LargePrimeFilterOptions {
            modulus: P,
            ..Default::default()
        };
        let (complete, report) = filter_large_primes(inputs.clone(), &options).unwrap();
        assert_eq!(complete.len(), 1);
        assert_eq!(report.cycles_closed, 1);
        let replay = replay_complete_relation(&complete[0], &inputs, P).unwrap();
        assert_eq!(replay.small, complete[0].small);
        assert_eq!(replay.rhs, complete[0].rhs);
    }

    #[test]
    fn double_large_prime_graph_cycle_closes() {
        let inputs = vec![
            row(1, &[(0, 1)], &[(10, 1), (11, P - 1)], 4),
            row(2, &[(1, 1)], &[(11, 1), (12, P - 1)], 7),
            row(3, &[(2, 1)], &[(12, 1), (10, P - 1)], 9),
        ];
        let options = LargePrimeFilterOptions {
            modulus: P,
            ..Default::default()
        };
        let (complete, report) = filter_large_primes(inputs.clone(), &options).unwrap();
        assert_eq!(complete.len(), 1, "{report:?}");
        assert_eq!(report.core_large_vertices, 3);
        assert!(report.merges >= 2);
        assert!(replay_complete_relation(&complete[0], &inputs, P).is_some());
    }

    #[test]
    fn hypergraph_peeling_removes_the_dangling_component() {
        let inputs = vec![
            row(1, &[(0, 1)], &[(10, 1), (11, 1), (12, 1)], 4),
            row(2, &[(1, 1)], &[(11, 1), (12, 1)], 7),
            // vertex 99 has degree one, so this row peels; that then peels
            // the first two as their remaining degrees collapse.
            row(3, &[(2, 1)], &[(12, 1), (99, 1)], 9),
        ];
        let options = LargePrimeFilterOptions {
            modulus: P,
            ..Default::default()
        };
        let (complete, report) = filter_large_primes(inputs, &options).unwrap();
        assert!(complete.is_empty());
        assert_eq!(report.peeled_relations, 3);
        assert_eq!(report.core_relations, 0);
    }

    #[test]
    fn arbitrary_coefficients_cancel_modulo_the_group_order() {
        let inputs = vec![
            row(1, &[(0, 3)], &[(42, 7)], 5),
            row(2, &[(1, 4)], &[(42, 9)], 6),
        ];
        let options = LargePrimeFilterOptions {
            modulus: P,
            ..Default::default()
        };
        let (complete, _) = filter_large_primes(inputs.clone(), &options).unwrap();
        assert_eq!(complete.len(), 1);
        let replay = replay_complete_relation(&complete[0], &inputs, P).unwrap();
        assert_eq!(replay, complete[0]);
    }

    #[test]
    fn a_fill_rejection_never_replaces_an_existing_pivot() {
        let inputs = vec![
            row(1, &[(0, 1)], &[(10, 1), (11, P - 1)], 4),
            row(2, &[(1, 1)], &[(11, 1), (12, P - 1)], 7),
            row(3, &[(2, 1)], &[(12, 1), (10, P - 1)], 9),
        ];
        let options = LargePrimeFilterOptions {
            modulus: P,
            max_small_weight: 1,
            ..Default::default()
        };
        let (complete, report) = filter_large_primes(inputs, &options).unwrap();
        assert!(complete.is_empty());
        assert_eq!(report.rejected_fill, 1);
        assert_eq!(
            report.pivots, 2,
            "the rejected row must not replace a pivot"
        );
    }
}
