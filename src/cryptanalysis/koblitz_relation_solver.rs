//! # Incremental modular echelon form for index-calculus relation matrices.
//!
//! The Koblitz index-calculus driver collects relations
//!
//! ```text
//!     Σ_o c_o · x_o  −  (h·b)·d  ≡  h·a   (mod r)
//! ```
//!
//! one at a time and wants to know, as early as possible, whether the
//! target scalar `d` is already pinned down.  Re-running a dense
//! `BigUint` Gaussian elimination after every batch is `O(U³)` big-integer
//! work per check and, worse, cannot distinguish "not yet determined"
//! from "determined": it returns a partial solution that must then be
//! verified in the group.
//!
//! This solver keeps the relation matrix in **reduced row echelon form
//! over `Z/rZ`** and updates it per row in `O(U²)` native `u64`
//! operations when `r < 2^64`. The wide subgroup path uses the same
//! echelon invariant with full-width coefficients. Column `U` is `d`;
//! because every pivot row has zeros in all other pivot columns, a pivot
//! in column `U` reads off `d`
//! directly, so `target()` answers the "is `d` determined?" question
//! exactly.  Dependent relations are recognised and dropped instead of
//! silently padding the matrix, and an inconsistent row is reported —
//! it can only come from a wrong relation, which the driver treats as a
//! failure rather than something to average away.
//!
//! The dense solver in
//! [`crate::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n`]
//! remains the reference; the tests cross-check the two on random
//! systems with planted solutions.

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::koblitz_index_calculus::KoblitzRelation;
use crate::utils::mod_inverse;

/// Outcome of feeding one relation to the solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RowStatus {
    /// The row raised the rank.
    Independent,
    /// The row was a linear combination of earlier rows; nothing changed.
    Dependent,
    /// The row contradicts earlier rows: `0 ≡ c` with `c ≠ 0`.  Only a
    /// wrong relation can do this, so the driver must fail closed.
    Inconsistent,
}

/// Reduced row echelon form of the relation matrix, maintained one row
/// at a time.  Columns `0..unknowns` are orbit logarithms, column
/// `unknowns` is the target scalar `d`, and the right-hand side sits in
/// the last position of each stored row.
#[derive(Clone, Debug)]
pub struct IncrementalRelationSolver {
    modulus: u64,
    unknowns: usize,
    /// `pivot_rows[c]` is the row whose leading entry (normalised to 1)
    /// is in column `c`, or `None` if column `c` is free.
    pivot_rows: Vec<Option<Vec<u64>>>,
    rank: usize,
    rows_seen: usize,
    dependent_rows: usize,
    inconsistent: bool,
}

#[inline]
fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

#[inline]
fn submod(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + (m - b)
    }
}

/// Modular inverse by the extended Euclidean algorithm; `None` when
/// `a` and `m` are not coprime (impossible for prime `m` and `a ≠ 0`).
fn invmod(a: u64, m: u64) -> Option<u64> {
    let (mut old_r, mut r) = (a as i128, m as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    if old_r != 1 {
        return None;
    }
    Some(old_s.rem_euclid(m as i128) as u64)
}

/// Export for the factor-base-logs precompute, which tracks coefficient
/// rank in `u64` while its relation rows stay big-integer for the final
/// dense solve.
pub(crate) fn to_u64_mod(v: &BigUint, m: u64) -> u64 {
    let reduced = v % BigUint::from(m);
    reduced.to_u64_digits().first().copied().unwrap_or(0)
}

impl IncrementalRelationSolver {
    /// A solver over `Z/rZ` with `unknowns` orbit columns plus the
    /// target column.  Returns `None` when the modulus does not fit a
    /// `u64` or is not at least 2; callers then fall back to the dense
    /// big-integer path.
    pub fn new(unknowns: usize, modulus: &BigUint) -> Option<Self> {
        if modulus.bits() > 64 || modulus < &BigUint::from(2u32) {
            return None;
        }
        let modulus = modulus.to_u64_digits().first().copied()?;
        Some(Self {
            modulus,
            unknowns,
            pivot_rows: vec![None; unknowns + 1],
            rank: 0,
            rows_seen: 0,
            dependent_rows: 0,
            inconsistent: false,
        })
    }

    /// The prime modulus `r`.
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Current rank of the relation matrix (independent rows kept).
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Relations fed so far, dependent ones included.
    pub fn rows_seen(&self) -> usize {
        self.rows_seen
    }

    /// Relations that were linear combinations of earlier ones.
    pub fn dependent_rows(&self) -> usize {
        self.dependent_rows
    }

    /// Whether an inconsistent row has ever been seen.
    pub fn is_inconsistent(&self) -> bool {
        self.inconsistent
    }

    /// Number of orbit columns.
    pub fn unknowns(&self) -> usize {
        self.unknowns
    }

    /// Feed a relation `[a]G + [b]Q = Σ P_i`, already rewritten over the
    /// orbit columns, with the cofactor `h` applied as in the driver.
    pub fn add_relation(&mut self, relation: &KoblitzRelation, cofactor: &BigUint) -> RowStatus {
        let m = self.modulus;
        let h = to_u64_mod(cofactor, m);
        let hb = mulmod(h, to_u64_mod(&relation.coef_b, m), m);
        let ha = mulmod(h, to_u64_mod(&relation.coef_a, m), m);
        let mut row = vec![0u64; self.unknowns + 2];
        for (o, c) in relation.row.iter().enumerate().take(self.unknowns) {
            row[o] = to_u64_mod(c, m);
        }
        row[self.unknowns] = submod(0, hb, m);
        row[self.unknowns + 1] = ha;
        self.add_row(row)
    }

    /// Feed a raw row `[c_0, …, c_{U−1}, c_d | rhs]` reduced mod `r`.
    pub fn add_row(&mut self, mut row: Vec<u64>) -> RowStatus {
        assert_eq!(row.len(), self.unknowns + 2, "row width");
        let m = self.modulus;
        let width = self.unknowns + 2;
        self.rows_seen += 1;
        // Reduce against every existing pivot, in column order.  A pivot
        // row has zeros before its pivot column, so subtraction starts
        // there.
        for col in 0..=self.unknowns {
            let factor = row[col];
            if factor == 0 {
                continue;
            }
            if let Some(pivot) = &self.pivot_rows[col] {
                for k in col..width {
                    if pivot[k] != 0 {
                        row[k] = submod(row[k], mulmod(factor, pivot[k], m), m);
                    }
                }
            }
        }
        let lead = (0..=self.unknowns).find(|&c| row[c] != 0);
        let Some(lead) = lead else {
            if row[self.unknowns + 1] != 0 {
                self.inconsistent = true;
                return RowStatus::Inconsistent;
            }
            self.dependent_rows += 1;
            return RowStatus::Dependent;
        };
        // Normalise, then clear this column from every other pivot row
        // so the form stays fully reduced.
        let inv = invmod(row[lead], m).expect("prime modulus: nonzero entries invert");
        for k in lead..width {
            row[k] = mulmod(row[k], inv, m);
        }
        for col in 0..=self.unknowns {
            if col == lead {
                continue;
            }
            if let Some(pivot) = self.pivot_rows[col].as_mut() {
                let factor = pivot[lead];
                if factor == 0 {
                    continue;
                }
                for k in lead..width {
                    if row[k] != 0 {
                        pivot[k] = submod(pivot[k], mulmod(factor, row[k], m), m);
                    }
                }
            }
        }
        self.pivot_rows[lead] = Some(row);
        self.rank += 1;
        RowStatus::Independent
    }

    /// The target scalar `d`, exactly when the relations determine it —
    /// i.e. the reduced form has a pivot in the `d` column.
    pub fn target(&self) -> Option<u64> {
        let row = self.pivot_rows[self.unknowns].as_ref()?;
        // In reduced form the pivot row for the last unknown column is
        // `[0, …, 0, 1 | rhs]`.
        debug_assert!(row[..self.unknowns].iter().all(|&c| c == 0));
        Some(row[self.unknowns + 1])
    }

    /// The target scalar as a `BigUint`, if determined.
    pub fn target_biguint(&self) -> Option<BigUint> {
        self.target().map(BigUint::from)
    }
}

/// Full-width reduced echelon form for a prime subgroup order above `u64`.
/// A target is returned only when its own column has a pivot. This avoids
/// treating the dense reference solver's arbitrary free-column values as a
/// rank certificate.
#[derive(Clone, Debug)]
pub struct WideIncrementalRelationSolver {
    modulus: BigUint,
    unknowns: usize,
    pivot_rows: Vec<Option<Vec<BigUint>>>,
    rank: usize,
    rows_seen: usize,
    dependent_rows: usize,
    inconsistent: bool,
}

impl WideIncrementalRelationSolver {
    pub fn new(unknowns: usize, modulus: &BigUint) -> Option<Self> {
        if modulus < &BigUint::from(2u32) {
            return None;
        }
        Some(Self {
            modulus: modulus.clone(),
            unknowns,
            pivot_rows: vec![None; unknowns + 1],
            rank: 0,
            rows_seen: 0,
            dependent_rows: 0,
            inconsistent: false,
        })
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn rows_seen(&self) -> usize {
        self.rows_seen
    }

    pub fn dependent_rows(&self) -> usize {
        self.dependent_rows
    }

    pub fn is_inconsistent(&self) -> bool {
        self.inconsistent
    }

    pub fn add_relation(&mut self, relation: &KoblitzRelation, cofactor: &BigUint) -> RowStatus {
        assert_eq!(relation.row.len(), self.unknowns, "relation width");
        let m = &self.modulus;
        let h = cofactor % m;
        let mut row = relation.row.clone();
        row.push((m - (&h * &relation.coef_b) % m) % m);
        row.push((&h * &relation.coef_a) % m);
        self.add_row(row)
    }

    /// Feed `[c_0, …, c_{U−1}, c_d | rhs]`, reducing every input modulo `r`.
    pub fn add_row(&mut self, mut row: Vec<BigUint>) -> RowStatus {
        assert_eq!(row.len(), self.unknowns + 2, "row width");
        let m = &self.modulus;
        let width = row.len();
        for value in &mut row {
            *value %= m;
        }
        self.rows_seen += 1;
        for col in 0..=self.unknowns {
            if row[col].is_zero() {
                continue;
            }
            if let Some(pivot) = &self.pivot_rows[col] {
                let factor = row[col].clone();
                for k in col..width {
                    if !pivot[k].is_zero() {
                        let term = (&factor * &pivot[k]) % m;
                        row[k] = (&row[k] + m - term) % m;
                    }
                }
            }
        }
        let lead = (0..=self.unknowns).find(|&col| !row[col].is_zero());
        let Some(lead) = lead else {
            if !row[self.unknowns + 1].is_zero() {
                self.inconsistent = true;
                return RowStatus::Inconsistent;
            }
            self.dependent_rows += 1;
            return RowStatus::Dependent;
        };
        let inverse = mod_inverse(&row[lead], m).expect("prime subgroup order");
        for value in &mut row[lead..width] {
            *value = (&*value * &inverse) % m;
        }
        for col in 0..=self.unknowns {
            if col == lead {
                continue;
            }
            if let Some(pivot) = self.pivot_rows[col].as_mut() {
                let factor = pivot[lead].clone();
                if factor.is_zero() {
                    continue;
                }
                for k in lead..width {
                    if !row[k].is_zero() {
                        let term = (&factor * &row[k]) % m;
                        pivot[k] = (&pivot[k] + m - term) % m;
                    }
                }
            }
        }
        self.pivot_rows[lead] = Some(row);
        self.rank += 1;
        RowStatus::Independent
    }

    pub fn target_biguint(&self) -> Option<BigUint> {
        let row = self.pivot_rows[self.unknowns].as_ref()?;
        debug_assert!(row[..self.unknowns].iter().all(BigUint::is_zero));
        debug_assert_eq!(row[self.unknowns], BigUint::one());
        Some(row[self.unknowns + 1].clone())
    }
}

/// Select native-word arithmetic when possible and full-width arithmetic
/// otherwise, with one rank/target interface for the index-calculus driver.
#[derive(Clone, Debug)]
pub enum RelationSolver {
    Word(IncrementalRelationSolver),
    Wide(WideIncrementalRelationSolver),
}

impl RelationSolver {
    pub fn new(unknowns: usize, modulus: &BigUint) -> Option<Self> {
        if modulus.bits() <= 64 {
            IncrementalRelationSolver::new(unknowns, modulus).map(Self::Word)
        } else {
            WideIncrementalRelationSolver::new(unknowns, modulus).map(Self::Wide)
        }
    }

    pub fn rank(&self) -> usize {
        match self {
            Self::Word(solver) => solver.rank(),
            Self::Wide(solver) => solver.rank(),
        }
    }

    pub fn add_relation(&mut self, relation: &KoblitzRelation, cofactor: &BigUint) -> RowStatus {
        match self {
            Self::Word(solver) => solver.add_relation(relation, cofactor),
            Self::Wide(solver) => solver.add_relation(relation, cofactor),
        }
    }

    pub fn target_biguint(&self) -> Option<BigUint> {
        match self {
            Self::Word(solver) => solver.target_biguint(),
            Self::Wide(solver) => solver.target_biguint(),
        }
    }
}

/// Minimal rank tracker for dense `u64` rows over `Z/mZ`.
///
/// The factor-base-logs precompute grows a big-integer relation matrix
/// and re-runs a dense solve after (almost) every new row, but a solve
/// can only succeed once the coefficient rows reach full column rank.
/// Tracking that rank here in `O(cols²)` native operations per row lets
/// the driver skip every provably-doomed dense attempt; the calls it
/// still makes — and their results — are unchanged.
#[derive(Clone, Debug, Default)]
pub struct U64RankTracker {
    modulus: u64,
    cols: usize,
    /// Echelon rows with strictly increasing leading positions.
    basis: Vec<(usize, Vec<u64>)>,
}

impl U64RankTracker {
    /// Track rank over `Z/mZ` for rows of exactly `cols` entries.
    pub fn new(modulus: u64, cols: usize) -> Self {
        Self {
            modulus,
            cols,
            basis: Vec::new(),
        }
    }

    /// Current rank (independent rows kept).
    pub fn rank(&self) -> usize {
        self.basis.len()
    }

    /// Insert a row; returns the new rank.  Rows are reduced against the
    /// basis in leading-position order, and newly kept rows are
    /// normalised to a monic leading entry, so a surviving nonzero row
    /// always carries a fresh leading position and the basis stays
    /// sorted.  One modular inverse per independent row.
    pub fn insert(&mut self, mut row: Vec<u64>) -> usize {
        assert_eq!(row.len(), self.cols, "row width");
        let m = self.modulus;
        for (lead, brow) in &self.basis {
            let factor = row[*lead];
            if factor == 0 {
                continue;
            }
            for k in *lead..self.cols {
                if brow[k] != 0 {
                    row[k] = submod(row[k], mulmod(factor, brow[k], m), m);
                }
            }
        }
        if let Some(lead) = row.iter().position(|&c| c != 0) {
            let inv = invmod(row[lead], m).expect("prime modulus: nonzero entries invert");
            if inv != 1 {
                for k in lead..self.cols {
                    row[k] = mulmod(row[k], inv, m);
                }
            }
            let pos = self
                .basis
                .iter()
                .position(|(l, _)| *l > lead)
                .unwrap_or(self.basis.len());
            self.basis.insert(pos, (lead, row));
        }
        self.basis.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n;
    use num_traits::Zero;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    fn planted_system(
        rng: &mut StdRng,
        unknowns: usize,
        rows: usize,
        modulus: u64,
        density: usize,
    ) -> (Vec<Vec<u64>>, Vec<u64>) {
        let secret: Vec<u64> = (0..=unknowns).map(|_| rng.gen_range(0..modulus)).collect();
        let mut out = Vec::new();
        for _ in 0..rows {
            let mut row = vec![0u64; unknowns + 2];
            for _ in 0..density {
                let o = rng.gen_range(0..unknowns);
                row[o] = rng.gen_range(1..modulus);
            }
            row[unknowns] = rng.gen_range(1..modulus);
            let mut rhs = 0u64;
            for (k, &c) in row.iter().enumerate().take(unknowns + 1) {
                rhs = (rhs + mulmod(c, secret[k], modulus)) % modulus;
            }
            row[unknowns + 1] = rhs;
            out.push(row);
        }
        (out, secret)
    }

    #[test]
    fn determines_the_planted_target_and_agrees_with_the_dense_solver() {
        let mut rng = StdRng::seed_from_u64(0x1234);
        for modulus in [127u64, 991, 2003, 262_543, 4_294_967_291] {
            for unknowns in [1usize, 3, 7, 12] {
                let (rows, secret) = planted_system(&mut rng, unknowns, unknowns + 6, modulus, 2);
                let big = BigUint::from(modulus);
                let mut solver = IncrementalRelationSolver::new(unknowns, &big).unwrap();
                let mut determined_at = None;
                for (i, row) in rows.iter().enumerate() {
                    let status = solver.add_row(row.clone());
                    assert_ne!(status, RowStatus::Inconsistent);
                    if let Some(d) = solver.target() {
                        assert_eq!(d, secret[unknowns], "modulus {modulus}, U {unknowns}");
                        if determined_at.is_none() {
                            determined_at = Some(i + 1);
                        }
                    }
                    // Whenever the dense solver's answer verifies, ours
                    // must be determined too, and equal.
                    let mut matrix: Vec<Vec<BigUint>> = rows[..=i]
                        .iter()
                        .map(|r| r[..=unknowns].iter().map(|&c| BigUint::from(c)).collect())
                        .collect();
                    let mut rhs: Vec<BigUint> = rows[..=i]
                        .iter()
                        .map(|r| BigUint::from(r[unknowns + 1]))
                        .collect();
                    if let Some(sol) = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, &big) {
                        let dense_d = sol[unknowns].to_u64_digits().first().copied().unwrap_or(0);
                        if solver.target().is_some() {
                            assert_eq!(dense_d, solver.target().unwrap());
                        }
                    }
                }
                assert!(
                    determined_at.is_some(),
                    "U + 6 random rows should determine d (modulus {modulus}, U {unknowns})"
                );
                assert!(solver.rank() <= unknowns + 1);
                assert_eq!(solver.rows_seen(), unknowns + 6);
                assert_eq!(solver.rank() + solver.dependent_rows(), unknowns + 6);
            }
        }
    }

    #[test]
    fn dependent_rows_do_not_raise_rank_and_wrong_rows_are_flagged() {
        let modulus = 991u64;
        let big = BigUint::from(modulus);
        let mut solver = IncrementalRelationSolver::new(3, &big).unwrap();
        let row = vec![5u64, 0, 7, 3, 11];
        assert_eq!(solver.add_row(row.clone()), RowStatus::Independent);
        // Twice the same row is dependent.
        let doubled: Vec<u64> = row.iter().map(|&c| (2 * c) % modulus).collect();
        assert_eq!(solver.add_row(doubled), RowStatus::Dependent);
        assert_eq!(solver.rank(), 1);
        // Same coefficients, different right-hand side: inconsistent.
        let mut wrong = row.clone();
        wrong[4] = 12;
        assert_eq!(solver.add_row(wrong), RowStatus::Inconsistent);
        assert!(solver.is_inconsistent());
        assert!(solver.target().is_none());
    }

    #[test]
    fn target_is_only_reported_when_pinned() {
        // Two unknowns and d; one row cannot determine d.
        let modulus = 127u64;
        let big = BigUint::from(modulus);
        let mut solver = IncrementalRelationSolver::new(2, &big).unwrap();
        assert_eq!(solver.add_row(vec![1, 1, 1, 5]), RowStatus::Independent);
        assert!(solver.target().is_none());
        assert_eq!(solver.add_row(vec![1, 0, 1, 3]), RowStatus::Independent);
        assert!(solver.target().is_none());
        assert_eq!(solver.add_row(vec![0, 1, 1, 4]), RowStatus::Independent);
        // x0 + x1 + d = 5, x0 + d = 3, x1 + d = 4  ⇒  x0 = 1, x1 = 2, d = 2.
        assert_eq!(solver.target(), Some(2));
    }

    #[test]
    fn rejects_moduli_that_do_not_fit() {
        let too_big = BigUint::from(1u8) << 70;
        assert!(IncrementalRelationSolver::new(4, &too_big).is_none());
        assert!(IncrementalRelationSolver::new(4, &BigUint::zero()).is_none());
    }

    #[test]
    fn wide_rank_and_target_use_the_entire_subgroup_order() {
        let modulus = BigUint::parse_bytes(b"2417851639230796216685689", 10).unwrap();
        let cofactor = BigUint::from(4u32);
        let inverse_h = mod_inverse(&cofactor, &modulus).unwrap();
        let column_log = BigUint::from(1234u32);
        let target_log = BigUint::from(5678u32);
        let high = BigUint::one() << 75usize;
        let first = &high + BigUint::from(9u32);
        let second = &first * BigUint::from(2u32) + BigUint::one();
        let make_relation = |coefficient: BigUint, b: u32| {
            let b = BigUint::from(b);
            let rhs = (&coefficient * &column_log + &modulus
                - (&cofactor * &b * &target_log) % &modulus)
                % &modulus;
            KoblitzRelation {
                coef_a: (rhs * &inverse_h) % &modulus,
                coef_b: b,
                summands: Vec::new(),
                summand_negated: Vec::new(),
                row: vec![coefficient],
            }
        };
        let relation1 = make_relation(first, 2);
        let relation2 = make_relation(second, 5);
        let mut solver = RelationSolver::new(1, &modulus).unwrap();
        assert!(matches!(&solver, RelationSolver::Wide(_)));
        assert_eq!(
            solver.add_relation(&relation1, &cofactor),
            RowStatus::Independent
        );
        assert_eq!(solver.rank(), 1);
        assert!(solver.target_biguint().is_none());
        assert_eq!(
            solver.add_relation(&relation2, &cofactor),
            RowStatus::Independent
        );
        assert_eq!(solver.rank(), 2);
        assert_eq!(solver.target_biguint(), Some(target_log.clone()));
        assert_eq!(
            solver.add_relation(&relation1, &cofactor),
            RowStatus::Dependent
        );
        let mut wrong = relation1.clone();
        wrong.coef_a = (&wrong.coef_a + BigUint::one()) % &modulus;
        assert_eq!(
            solver.add_relation(&wrong, &cofactor),
            RowStatus::Inconsistent
        );
        let RelationSolver::Wide(wide) = solver else {
            unreachable!();
        };
        assert_eq!(wide.rows_seen(), 4);
        assert_eq!(wide.dependent_rows(), 1);
        assert!(wide.is_inconsistent());

        let mut matrix = vec![
            vec![
                relation1.row[0].clone(),
                (&modulus - (&cofactor * &relation1.coef_b) % &modulus) % &modulus,
            ],
            vec![
                relation2.row[0].clone(),
                (&modulus - (&cofactor * &relation2.coef_b) % &modulus) % &modulus,
            ],
        ];
        let mut rhs = vec![
            (&cofactor * &relation1.coef_a) % &modulus,
            (&cofactor * &relation2.coef_a) % &modulus,
        ];
        let dense = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, &modulus).unwrap();
        assert_eq!(dense, vec![column_log, target_log]);
    }

    #[test]
    fn wide_random_planted_rows_agree_with_dense_solution() {
        let modulus = BigUint::parse_bytes(b"2417851639230796216685689", 10).unwrap();
        let mut rng = StdRng::seed_from_u64(0x83_1c_2026);
        for unknowns in [1usize, 3, 6] {
            let columns = unknowns + 1;
            let mut sample = || {
                (BigUint::from(rng.gen::<u64>()) + (BigUint::from(rng.gen::<u64>()) << 64usize))
                    % &modulus
            };
            let secret: Vec<BigUint> = (0..columns).map(|_| sample()).collect();
            let mut solver = WideIncrementalRelationSolver::new(unknowns, &modulus).unwrap();
            let mut matrix = Vec::new();
            let mut rhs_values = Vec::new();
            for _ in 0..(2 * columns + 5) {
                let coefficients: Vec<BigUint> = (0..columns).map(|_| sample()).collect();
                let rhs = coefficients
                    .iter()
                    .zip(&secret)
                    .fold(BigUint::zero(), |acc, (coefficient, value)| {
                        (acc + coefficient * value) % &modulus
                    });
                let mut row = coefficients.clone();
                row.push(rhs.clone());
                assert_ne!(solver.add_row(row), RowStatus::Inconsistent);
                if let Some(target) = solver.target_biguint() {
                    assert_eq!(target, secret[unknowns]);
                }
                matrix.push(coefficients);
                rhs_values.push(rhs);
            }
            assert_eq!(solver.rank(), columns);
            assert_eq!(solver.target_biguint(), Some(secret[unknowns].clone()));
            assert_eq!(
                gaussian_eliminate_mod_n(&mut matrix, &mut rhs_values, &modulus).unwrap(),
                secret
            );
        }
    }

    #[test]
    fn rank_tracker_agrees_with_incremental_solver() {
        use rand::{rngs::StdRng, Rng, SeedableRng};
        // Same row stream into both structures: ranks must agree at
        // every step (the incremental solver also folds in a `d` column
        // and rhs, which can only keep its rank at or above the
        // coefficient rank — so assert exactly that).
        let m = 2003u64;
        let big = BigUint::from(m);
        let mut rng = StdRng::seed_from_u64(0x274c);
        for unknowns in [1usize, 4, 9] {
            let mut tracker = U64RankTracker::new(m, unknowns);
            let mut solver = IncrementalRelationSolver::new(unknowns, &big).expect("modulus fits");
            for _ in 0..3 * unknowns + 5 {
                let coeff: Vec<u64> = (0..unknowns).map(|_| rng.gen_range(0..m)).collect();
                let t_rank = tracker.insert(coeff.clone());
                let mut full = coeff;
                full.push(rng.gen_range(0..m));
                full.push(rng.gen_range(0..m));
                solver.add_row(full);
                assert!(solver.rank() >= t_rank, "U={unknowns}");
                assert!(t_rank <= unknowns, "U={unknowns}");
            }
            assert_eq!(tracker.rank(), unknowns, "U={unknowns} reaches full rank");
        }
    }

    #[test]
    fn rank_tracker_counts_basis_size() {
        use rand::{rngs::StdRng, Rng, SeedableRng};
        let m = 2003u64;
        // Identity rows: rank grows one per row.
        let mut tracker = U64RankTracker::new(m, 4);
        for i in 0..4 {
            let mut row = vec![0u64; 4];
            row[i] = 1;
            assert_eq!(tracker.insert(row), i + 1);
        }
        // Duplicates, multiples, and zero rows add nothing.
        assert_eq!(tracker.insert(vec![1, 0, 0, 0]), 4);
        assert_eq!(tracker.insert(vec![0, 0, 5, 0]), 4);
        assert_eq!(tracker.insert(vec![0, 0, 0, 0]), 4);
        // Full-rank random matrices reach min(rows, cols); extra rows stop.
        let mut rng = StdRng::seed_from_u64(0x7a);
        for cols in [1usize, 3, 7] {
            let mut tracker = U64RankTracker::new(m, cols);
            for i in 0..2 * cols + 3 {
                let row: Vec<u64> = (0..cols).map(|_| rng.gen_range(0..m)).collect();
                let rank = tracker.insert(row);
                assert!(rank <= cols.min(i + 1));
            }
            assert_eq!(tracker.rank(), cols, "random rows reach full rank");
        }
    }
}
