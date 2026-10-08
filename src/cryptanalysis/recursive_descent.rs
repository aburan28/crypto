//! Budgeted recursive individual-log descent.
//!
//! The scheduler is deliberately independent of a curve representation.  An
//! oracle supplies certified relations
//!
//! ```text
//! a * log(target) = constant + sum b_i * log(child_i)  (mod r),
//! ```
//!
//! and a canonicalisation map.  Canonicalisation is where signed Frobenius
//! orbits live: if `P = ±pi^k(P0)` and Frobenius acts as `lambda`, return
//! `(P0, ±lambda^k)`.  Ordinary recursive point decomposition returns the
//! identity multiplier, while large-prime descent simply returns large-prime
//! points as children.  Every relation is verified by the oracle before the
//! scheduler consumes it, and every final logarithm is independently verified.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;

#[inline]
fn add_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 + b as u128) % modulus as u128) as u64
}

#[inline]
fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
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

/// A canonical representative and the scalar transporting its logarithm back
/// to the original node: `log(original) = multiplier * log(node)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalNode<N> {
    pub node: N,
    pub multiplier: u64,
}

/// One certified decomposition candidate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DescentRelation<N, W> {
    /// Coefficient of `log(target)` on the left side.
    pub target_coefficient: u64,
    /// Known logarithmic constant on the right side.
    pub constant: u64,
    pub children: Vec<(N, u64)>,
    /// Collector-specific witness retained in the transcript.
    pub witness: W,
}

/// Curve- and oracle-specific operations used by [`solve_recursive_log`].
pub trait RecursiveDescentOracle {
    type Node: Clone + Eq + Hash;
    type Witness: Clone;

    /// Group order of the prime-order subgroup.
    fn modulus(&self) -> u64;

    /// Quotient a point by any proven symmetry.  Use multiplier one when no
    /// quotient applies.
    fn canonicalise(&self, node: &Self::Node) -> CanonicalNode<Self::Node>;

    /// Return a precomputed factor-base logarithm, if this is a leaf.
    fn base_log(&self, node: &Self::Node) -> Option<u64>;

    /// Produce randomized/recursive decomposition candidates.  `attempt`
    /// changes the oracle's deterministic seed and makes randomized smoothing
    /// replayable.
    fn relations(
        &mut self,
        node: &Self::Node,
        attempt: usize,
    ) -> Vec<DescentRelation<Self::Node, Self::Witness>>;

    /// Re-add the points (or replay the algebraic witness) before accepting a
    /// relation.  The scheduler never trusts a solver result without this.
    fn verify_relation(
        &self,
        node: &Self::Node,
        relation: &DescentRelation<Self::Node, Self::Witness>,
    ) -> bool;

    /// Independent final check, normally `[log]G == node`.
    fn verify_log(&self, node: &Self::Node, logarithm: u64) -> bool;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RecursiveDescentOptions {
    pub max_depth: usize,
    pub max_nodes: usize,
    pub attempts_per_node: usize,
    pub max_relations: usize,
    pub max_children_per_relation: usize,
}

impl Default for RecursiveDescentOptions {
    fn default() -> Self {
        Self {
            max_depth: 32,
            max_nodes: 100_000,
            attempts_per_node: 64,
            max_relations: 1_000_000,
            max_children_per_relation: 64,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecursiveDescentReport {
    pub nodes_visited: usize,
    pub base_hits: usize,
    pub memo_hits: usize,
    pub orbit_rewrites: usize,
    pub relations_tried: usize,
    pub invalid_relations: usize,
    pub noninvertible_relations: usize,
    pub cyclic_relations: usize,
    pub failed_children: usize,
    pub max_depth_reached: usize,
    pub final_verified: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RecursiveDescentFailure {
    InvalidModulus,
    NodeBudget,
    RelationBudget,
    DepthBudget,
    NoDecomposition,
    FinalVerification,
}

/// One accepted step in the replay transcript.
#[derive(Clone, Debug)]
pub struct DescentStep<N, W> {
    pub target: N,
    pub relation: DescentRelation<N, W>,
    pub child_logs: Vec<u64>,
    pub result: u64,
}

#[derive(Clone, Debug)]
pub struct RecursiveDescentSolution<N, W> {
    pub logarithm: u64,
    pub transcript: Vec<DescentStep<N, W>>,
    pub report: RecursiveDescentReport,
}

struct Scheduler<'a, O: RecursiveDescentOracle> {
    oracle: &'a mut O,
    options: RecursiveDescentOptions,
    modulus: u64,
    memo: HashMap<O::Node, u64>,
    active: HashSet<O::Node>,
    transcript: Vec<DescentStep<O::Node, O::Witness>>,
    report: RecursiveDescentReport,
    budget_failure: Option<RecursiveDescentFailure>,
}

impl<O: RecursiveDescentOracle> Scheduler<'_, O> {
    fn descend(&mut self, original: O::Node, depth: usize) -> Option<u64> {
        self.report.max_depth_reached = self.report.max_depth_reached.max(depth);
        if depth > self.options.max_depth {
            self.budget_failure
                .get_or_insert(RecursiveDescentFailure::DepthBudget);
            return None;
        }
        let canonical = self.oracle.canonicalise(&original);
        let multiplier = canonical.multiplier % self.modulus;
        if multiplier != 1 || canonical.node != original {
            self.report.orbit_rewrites += 1;
        }
        let node = canonical.node;
        if let Some(&logarithm) = self.memo.get(&node) {
            self.report.memo_hits += 1;
            return Some(mul_mod(multiplier, logarithm, self.modulus));
        }
        if self.report.nodes_visited >= self.options.max_nodes {
            self.budget_failure
                .get_or_insert(RecursiveDescentFailure::NodeBudget);
            return None;
        }
        self.report.nodes_visited += 1;
        if let Some(logarithm) = self.oracle.base_log(&node) {
            let logarithm = logarithm % self.modulus;
            if self.oracle.verify_log(&node, logarithm) {
                self.report.base_hits += 1;
                self.memo.insert(node, logarithm);
                return Some(mul_mod(multiplier, logarithm, self.modulus));
            }
        }
        if !self.active.insert(node.clone()) {
            self.report.cyclic_relations += 1;
            return None;
        }

        let mut solved = None;
        'attempts: for attempt in 0..self.options.attempts_per_node.max(1) {
            let relations = self.oracle.relations(&node, attempt);
            for relation in relations {
                if self.report.relations_tried >= self.options.max_relations {
                    self.budget_failure
                        .get_or_insert(RecursiveDescentFailure::RelationBudget);
                    break 'attempts;
                }
                self.report.relations_tried += 1;
                if relation.children.len() > self.options.max_children_per_relation
                    || !self.oracle.verify_relation(&node, &relation)
                {
                    self.report.invalid_relations += 1;
                    continue;
                }
                let Some(inv_target) =
                    inverse_mod(relation.target_coefficient % self.modulus, self.modulus)
                else {
                    self.report.noninvertible_relations += 1;
                    continue;
                };
                let mut rhs = relation.constant % self.modulus;
                let mut child_logs = Vec::with_capacity(relation.children.len());
                let mut all_children = true;
                for (child, coefficient) in &relation.children {
                    if self.active.contains(&self.oracle.canonicalise(child).node) {
                        self.report.cyclic_relations += 1;
                        all_children = false;
                        break;
                    }
                    let Some(logarithm) = self.descend(child.clone(), depth + 1) else {
                        self.report.failed_children += 1;
                        all_children = false;
                        break;
                    };
                    child_logs.push(logarithm);
                    rhs = add_mod(
                        rhs,
                        mul_mod(*coefficient % self.modulus, logarithm, self.modulus),
                        self.modulus,
                    );
                }
                if !all_children {
                    continue;
                }
                let logarithm = mul_mod(rhs, inv_target, self.modulus);
                if !self.oracle.verify_log(&node, logarithm) {
                    self.report.invalid_relations += 1;
                    continue;
                }
                self.transcript.push(DescentStep {
                    target: node.clone(),
                    relation,
                    child_logs,
                    result: logarithm,
                });
                solved = Some(logarithm);
                break 'attempts;
            }
        }
        self.active.remove(&node);
        if let Some(logarithm) = solved {
            self.memo.insert(node, logarithm);
            Some(mul_mod(multiplier, logarithm, self.modulus))
        } else {
            None
        }
    }
}

/// Solve one target logarithm by randomized recursive decomposition.
pub fn solve_recursive_log<O: RecursiveDescentOracle>(
    oracle: &mut O,
    target: O::Node,
    options: &RecursiveDescentOptions,
) -> Result<RecursiveDescentSolution<O::Node, O::Witness>, RecursiveDescentFailure> {
    let modulus = oracle.modulus();
    // Signed extended Euclid is bounded for the subgroup orders supported by
    // the ECDLP sparse stack.
    if !(2..(1u64 << 63)).contains(&modulus) {
        return Err(RecursiveDescentFailure::InvalidModulus);
    }
    let original = target.clone();
    let mut scheduler = Scheduler {
        oracle,
        options: *options,
        modulus,
        memo: HashMap::new(),
        active: HashSet::new(),
        transcript: Vec::new(),
        report: RecursiveDescentReport::default(),
        budget_failure: None,
    };
    let logarithm = scheduler.descend(target, 0).ok_or_else(|| {
        scheduler
            .budget_failure
            .clone()
            .unwrap_or(RecursiveDescentFailure::NoDecomposition)
    })?;
    if !scheduler.oracle.verify_log(&original, logarithm) {
        return Err(RecursiveDescentFailure::FinalVerification);
    }
    scheduler.report.final_verified = true;
    Ok(RecursiveDescentSolution {
        logarithm,
        transcript: scheduler.transcript,
        report: scheduler.report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct ToyOracle {
        logs: HashMap<u32, u64>,
        base: HashSet<u32>,
    }

    impl RecursiveDescentOracle for ToyOracle {
        type Node = u32;
        type Witness = &'static str;

        fn modulus(&self) -> u64 {
            101
        }

        fn canonicalise(&self, node: &u32) -> CanonicalNode<u32> {
            // Toy signed orbit: 100+n has log -log(n).
            if *node >= 100 {
                CanonicalNode {
                    node: *node - 100,
                    multiplier: 100,
                }
            } else {
                CanonicalNode {
                    node: *node,
                    multiplier: 1,
                }
            }
        }

        fn base_log(&self, node: &u32) -> Option<u64> {
            self.base.contains(node).then(|| self.logs[node])
        }

        fn relations(
            &mut self,
            node: &u32,
            _attempt: usize,
        ) -> Vec<DescentRelation<u32, &'static str>> {
            match *node {
                4 => vec![DescentRelation {
                    target_coefficient: 1,
                    constant: 3,
                    children: vec![(3, 2), (102, 1)],
                    witness: "large-prime split",
                }],
                3 => vec![DescentRelation {
                    target_coefficient: 2,
                    constant: 5,
                    children: vec![(1, 4), (2, 3)],
                    witness: "recursive split",
                }],
                _ => Vec::new(),
            }
        }

        fn verify_relation(
            &self,
            node: &u32,
            relation: &DescentRelation<u32, &'static str>,
        ) -> bool {
            let lhs = mul_mod(relation.target_coefficient, self.logs[node], 101);
            let mut rhs = relation.constant;
            for (child, coefficient) in &relation.children {
                let canonical = self.canonicalise(child);
                let child_log = mul_mod(canonical.multiplier, self.logs[&canonical.node], 101);
                rhs = add_mod(rhs, mul_mod(*coefficient, child_log, 101), 101);
            }
            lhs == rhs
        }

        fn verify_log(&self, node: &u32, logarithm: u64) -> bool {
            let canonical = self.canonicalise(node);
            logarithm == mul_mod(canonical.multiplier, self.logs[&canonical.node], 101)
        }
    }

    #[test]
    fn recursive_large_prime_and_signed_orbit_descent_verifies() {
        // Choose the logs to satisfy both relations:
        // 2*l3 = 5 + 4*l1 + 3*l2; l4 = 3 + 2*l3 - l2.
        let l1 = 7;
        let l2 = 11;
        let l3 = mul_mod(5 + 4 * l1 + 3 * l2, inverse_mod(2, 101).unwrap(), 101);
        let l4 = (3 + 2 * l3 + 101 - l2) % 101;
        let mut oracle = ToyOracle {
            logs: [(1, l1), (2, l2), (3, l3), (4, l4)].into_iter().collect(),
            base: [1, 2].into_iter().collect(),
        };
        let solution =
            solve_recursive_log(&mut oracle, 4, &RecursiveDescentOptions::default()).unwrap();
        assert_eq!(solution.logarithm, l4);
        assert_eq!(solution.transcript.len(), 2);
        assert!(solution.report.orbit_rewrites >= 1);
        assert!(solution.report.final_verified);
    }

    #[test]
    fn cyclic_relation_is_rejected_without_recursing_forever() {
        #[derive(Clone)]
        struct Cycle;
        impl RecursiveDescentOracle for Cycle {
            type Node = u8;
            type Witness = ();
            fn modulus(&self) -> u64 {
                101
            }
            fn canonicalise(&self, node: &u8) -> CanonicalNode<u8> {
                CanonicalNode {
                    node: *node,
                    multiplier: 1,
                }
            }
            fn base_log(&self, _: &u8) -> Option<u64> {
                None
            }
            fn relations(&mut self, node: &u8, _: usize) -> Vec<DescentRelation<u8, ()>> {
                vec![DescentRelation {
                    target_coefficient: 1,
                    constant: 0,
                    children: vec![(*node, 1)],
                    witness: (),
                }]
            }
            fn verify_relation(&self, _: &u8, _: &DescentRelation<u8, ()>) -> bool {
                true
            }
            fn verify_log(&self, _: &u8, _: u64) -> bool {
                false
            }
        }
        let mut oracle = Cycle;
        let options = RecursiveDescentOptions {
            attempts_per_node: 1,
            ..Default::default()
        };
        assert_eq!(
            solve_recursive_log(&mut oracle, 1, &options).unwrap_err(),
            RecursiveDescentFailure::NoDecomposition
        );
    }
}
