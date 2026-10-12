//! Pair accounting for the CDCL shell: which pairs of the pair space each
//! learned clause excludes, and how many of them were never sent to the leaf.
//!
//! The pair space is `{(a, b) : a ≤ b, a, b ∈ S}` for the leaf's candidate set
//! `S` (`m = |S|` abscissae).  A *pair-only* clause (every literal on a bit of
//! `X₁` or `X₂`) excludes the pairs whose bits falsify all its literals: its
//! **cube**.  Three numbers per clause:
//!
//! * `cube`: valid pairs in the cube, whatever their history;
//! * `marginal`: pairs of the cube not excluded by an earlier clause (the pair
//!   whose refutation produced the clause is one of them);
//! * `saved`: the part of `marginal` whose pair the leaf had **not** evaluated
//!   when the clause was learned, i.e. the pairs the clause removed without a
//!   leaf call.
//!
//! For a refuted target every pair is either evaluated or saved, so the sum of
//! `saved` over the learned clauses is at most `pairs − leaf visits`.

use crate::shell::Layout;
use crypto_lib::cryptanalysis::sat::Lit;
use std::collections::BTreeMap;

pub struct PairAcct {
    lay: Layout,
    m: usize,
    space: Vec<u64>,
    idx_of: Vec<u32>,
    evaluated: Vec<u64>,
    excluded: Vec<u64>,
    pub total_pairs: u64,
    pub evaluated_pairs: u64,
    pub excluded_pairs: u64,
    pub learned: u64,
    pub learned_pair_only: u64,
    pub learned_other: u64,
    pub learned_units: u64,
    pub cube_sum: u64,
    pub marginal_sum: u64,
    pub saved_sum: u64,
    pub max_cube: u64,
    pub max_saved: u64,
    /// Pair-only learned clauses with fewer than `2ℓ` literals.
    pub shorter_than_pair: u64,
    pub len_hist: BTreeMap<usize, u64>,
    /// Plain encoding: pairs excluded by the theory's own clause.
    pub theory_excluded: u64,
}

impl PairAcct {
    pub fn new(lay: Layout, space: &[u64]) -> Self {
        let m = space.len();
        let mut idx_of = vec![u32::MAX; 1usize << lay.l];
        for (i, &x) in space.iter().enumerate() {
            idx_of[x as usize] = i as u32;
        }
        let words = (m * m).div_ceil(64);
        Self {
            lay,
            m,
            space: space.to_vec(),
            idx_of,
            evaluated: vec![0; words],
            excluded: vec![0; words],
            total_pairs: (m as u64) * (m as u64 + 1) / 2,
            evaluated_pairs: 0,
            excluded_pairs: 0,
            learned: 0,
            learned_pair_only: 0,
            learned_other: 0,
            learned_units: 0,
            cube_sum: 0,
            marginal_sum: 0,
            saved_sum: 0,
            max_cube: 0,
            max_saved: 0,
            shorter_than_pair: 0,
            len_hist: BTreeMap::new(),
            theory_excluded: 0,
        }
    }

    fn pid(&self, a: u64, b: u64) -> usize {
        let (ia, ib) = (self.idx_of[a as usize], self.idx_of[b as usize]);
        debug_assert!(ia != u32::MAX && ib != u32::MAX && ia <= ib);
        ia as usize * self.m + ib as usize
    }

    fn test(v: &[u64], p: usize) -> bool {
        v[p / 64] >> (p % 64) & 1 == 1
    }

    fn set(v: &mut [u64], p: usize) {
        v[p / 64] |= 1u64 << (p % 64);
    }

    /// The leaf evaluated the pair `(a, b)`.
    pub fn on_eval(&mut self, a: u64, b: u64) {
        let p = self.pid(a, b);
        if !Self::test(&self.evaluated, p) {
            Self::set(&mut self.evaluated, p);
            self.evaluated_pairs += 1;
        }
    }

    /// The theory's own clause excluded the pair (encoding without conflict analysis).
    pub fn on_theory_exclusion(&mut self, a: u64, b: u64) {
        let p = self.pid(a, b);
        if !Self::test(&self.excluded, p) {
            Self::set(&mut self.excluded, p);
            self.excluded_pairs += 1;
            self.theory_excluded += 1;
        }
    }

    /// Candidates of one summand matching the fixed bits: indices into the space.
    fn matches(&self, mask: u64, val: u64) -> Vec<usize> {
        if mask == (1u64 << self.lay.l) - 1 {
            let i = self.idx_of[val as usize];
            return if i == u32::MAX { vec![] } else { vec![i as usize] };
        }
        self.space
            .iter()
            .enumerate()
            .filter(|(_, &x)| x & mask == val)
            .map(|(i, _)| i)
            .collect()
    }

    /// A clause conflict analysis learned.
    pub fn on_learnt(&mut self, clause: &[Lit]) {
        self.learned += 1;
        *self.len_hist.entry(clause.len()).or_insert(0) += 1;
        if clause.len() == 1 {
            self.learned_units += 1;
        }
        let l = self.lay.l;
        let (mut m1, mut v1, mut m2, mut v2) = (0u64, 0u64, 0u64, 0u64);
        for &lit in clause {
            let var = lit.unsigned_abs();
            // A literal is false when the variable takes the value that falsifies it.
            let falsifying_one = lit < 0;
            if (1..=l).contains(&var) {
                let j = l - 1 - (var - 1);
                m1 |= 1 << j;
                if falsifying_one {
                    v1 |= 1 << j;
                }
            } else if (l + 1..=2 * l).contains(&var) {
                let j = l - 1 - (var - l - 1);
                m2 |= 1 << j;
                if falsifying_one {
                    v2 |= 1 << j;
                }
            } else {
                self.learned_other += 1;
                return;
            }
        }
        self.learned_pair_only += 1;
        if clause.len() < 2 * l as usize {
            self.shorter_than_pair += 1;
        }
        let (l1, l2) = (self.matches(m1, v1), self.matches(m2, v2));
        let (mut cube, mut marginal, mut saved) = (0u64, 0u64, 0u64);
        for &i1 in &l1 {
            for &i2 in &l2 {
                if i1 > i2 {
                    continue;
                }
                cube += 1;
                let p = i1 * self.m + i2;
                if !Self::test(&self.excluded, p) {
                    Self::set(&mut self.excluded, p);
                    self.excluded_pairs += 1;
                    marginal += 1;
                    if !Self::test(&self.evaluated, p) {
                        saved += 1;
                    }
                }
            }
        }
        self.cube_sum += cube;
        self.marginal_sum += marginal;
        self.saved_sum += saved;
        self.max_cube = self.max_cube.max(cube);
        self.max_saved = self.max_saved.max(saved);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube_marginal_and_saved_of_a_short_clause() {
        let lay = Layout::new(3);
        // Space: all of V = 0..8.
        let space: Vec<u64> = (0..8).collect();
        let mut a = PairAcct::new(lay, &space);
        assert_eq!(a.total_pairs, 36);
        // Evaluate (0, 0) and (0, 1): both in the cube of "x1 = 0, x2 high bit 0".
        a.on_eval(0, 0);
        a.on_eval(0, 1);
        // Clause forbidding x1 = 000 and the top bit of x2 = 0: literals on the variables of
        // x1 (vars 1..=3, msb first) and x2's msb (var 4): all false means x1 = 0, x2 msb = 0.
        let clause = [1, 2, 3, 4];
        a.on_learnt(&clause);
        // x1 = 0; x2 in {0,1,2,3} with x1 <= x2: four pairs; two evaluated.
        assert_eq!((a.cube_sum, a.marginal_sum, a.saved_sum), (4, 4, 2));
        // The same clause again adds a cube but no marginal pairs.
        a.on_learnt(&clause);
        assert_eq!((a.cube_sum, a.marginal_sum, a.saved_sum), (8, 4, 2));
        // A full-length clause for (3, 5): x1 = 011 (lits: var1 pos, var2 neg, var3 neg ... falsified by 0,1,1).
        let (x1, x2) = (3u64, 5u64);
        let mut c = Vec::new();
        for s in 0..2u32 {
            let x = if s == 0 { x1 } else { x2 };
            for j in (0..3).rev() {
                let var = (1 + s * 3 + (2 - j)) as i32;
                c.push(if x >> j & 1 == 1 { -var } else { var });
            }
        }
        a.on_learnt(&c);
        assert_eq!(a.learned_pair_only, 3);
        assert_eq!(a.marginal_sum, 5);
        assert_eq!(a.saved_sum, 3);
    }
}
