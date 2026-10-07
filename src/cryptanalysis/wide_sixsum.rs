//! Fixed-width Boolean system for the direct n83 six-summand S3 chain.
//! This admission path is independent of the 128-variable wide solver.

use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::fx_hash::FxMap;
use crate::cryptanalysis::gf2_elim;
use crate::cryptanalysis::wide_groebner::WideFieldTable;

pub const MAX_VARS: usize = 512;
pub const MAX_ROOT_COLS: usize = 1_500_000;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Mono512(pub [u64; 8]);

impl Mono512 {
    pub fn var(v: usize) -> Self {
        assert!(v < MAX_VARS);
        let mut words = [0; 8];
        words[v / 64] = 1u64 << (v % 64);
        Self(words)
    }

    pub fn degree(self) -> u32 {
        self.0.iter().map(|w| w.count_ones()).sum()
    }

    pub fn union(self, other: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] | other.0[i]))
    }

    pub fn divides_assignment(self, values: &Self) -> bool {
        self.0.iter().zip(values.0).all(|(a, b)| a & !b == 0)
    }

    fn intersects(self, other: Self) -> bool {
        self.0.iter().zip(other.0).any(|(a, b)| a & b != 0)
    }

    fn without(self, other: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] & !other.0[i]))
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Poly512 {
    pub terms: Vec<Mono512>,
}

impl Poly512 {
    pub fn zero() -> Self {
        Self::default()
    }

    pub fn one() -> Self {
        Self {
            terms: vec![Mono512::default()],
        }
    }

    pub fn var(v: usize) -> Self {
        Self {
            terms: vec![Mono512::var(v)],
        }
    }

    pub fn from_monos(mut terms: Vec<Mono512>) -> Self {
        terms.sort_unstable();
        let mut out = Vec::with_capacity(terms.len());
        let mut i = 0;
        while i < terms.len() {
            let mut j = i + 1;
            while j < terms.len() && terms[j] == terms[i] {
                j += 1;
            }
            if (j - i) & 1 == 1 {
                out.push(terms[i]);
            }
            i = j;
        }
        Self { terms: out }
    }

    pub fn add(&self, rhs: &Self) -> Self {
        let mut terms = Vec::with_capacity(self.terms.len() + rhs.terms.len());
        terms.extend_from_slice(&self.terms);
        terms.extend_from_slice(&rhs.terms);
        Self::from_monos(terms)
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut terms = Vec::with_capacity(self.terms.len().saturating_mul(rhs.terms.len()));
        for &a in &self.terms {
            for &b in &rhs.terms {
                terms.push(a.union(b));
            }
        }
        Self::from_monos(terms)
    }

    pub fn eval(&self, values: &Mono512) -> bool {
        self.terms
            .iter()
            .filter(|m| m.divides_assignment(values))
            .count()
            & 1
            == 1
    }

    pub fn degree(&self) -> u32 {
        self.terms.iter().map(|m| m.degree()).max().unwrap_or(0)
    }

    /// Substitute all selected Boolean variables in one pass. Terms
    /// containing a zero vanish; stripping one-bits can create pairs
    /// that must cancel in the Boolean quotient.
    pub fn assign_constants(&self, zeros: Mono512, ones: Mono512) -> Self {
        assert!(!zeros.intersects(ones));
        Self::from_monos(
            self.terms
                .iter()
                .copied()
                .filter(|m| !m.intersects(zeros))
                .map(|m| m.without(ones))
                .collect(),
        )
    }
}

#[derive(Clone)]
struct Sym512 {
    coords: Vec<Poly512>,
}

impl Sym512 {
    fn constant(value: &F2mElement, n: usize) -> Self {
        let words = value.raw_bits();
        Self {
            coords: (0..n)
                .map(|k| {
                    if words[k / 64] >> (k % 64) & 1 == 1 {
                        Poly512::one()
                    } else {
                        Poly512::zero()
                    }
                })
                .collect(),
        }
    }

    fn subspace(basis: &[F2mElement], offset: usize, n: usize) -> Self {
        let mut coords = vec![Vec::new(); n];
        for (i, base) in basis.iter().enumerate() {
            let words = base.raw_bits();
            for (k, terms) in coords.iter_mut().enumerate() {
                if words[k / 64] >> (k % 64) & 1 == 1 {
                    terms.push(Mono512::var(offset + i));
                }
            }
        }
        Self {
            coords: coords.into_iter().map(Poly512::from_monos).collect(),
        }
    }

    fn free(offset: usize, n: usize) -> Self {
        Self {
            coords: (0..n).map(|k| Poly512::var(offset + k)).collect(),
        }
    }

    fn add(&self, rhs: &Self) -> Self {
        Self {
            coords: self
                .coords
                .iter()
                .zip(&rhs.coords)
                .map(|(a, b)| a.add(b))
                .collect(),
        }
    }

    fn mul(&self, rhs: &Self, table: &impl WideFieldTable) -> Self {
        let n = table.degree() as usize;
        let mut acc: Vec<Vec<Mono512>> = vec![Vec::new(); n];
        for (i, a) in self.coords.iter().enumerate() {
            if a.terms.is_empty() {
                continue;
            }
            for (j, b) in rhs.coords.iter().enumerate() {
                if b.terms.is_empty() {
                    continue;
                }
                let product = a.mul(b);
                let mut bits = table.product_bits(i, j);
                while bits != 0 {
                    let k = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    acc[k].extend_from_slice(&product.terms);
                }
            }
        }
        Self {
            coords: acc.into_iter().map(Poly512::from_monos).collect(),
        }
    }

    fn square(&self, table: &impl WideFieldTable) -> Self {
        let n = table.degree() as usize;
        let mut acc: Vec<Vec<Mono512>> = vec![Vec::new(); n];
        for (i, a) in self.coords.iter().enumerate() {
            let mut bits = table.square_bits(i);
            while bits != 0 {
                let k = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                acc[k].extend_from_slice(&a.terms);
            }
        }
        Self {
            coords: acc.into_iter().map(Poly512::from_monos).collect(),
        }
    }
}

fn s3(
    x: &Sym512,
    y: &Sym512,
    z: &Sym512,
    b: &F2mElement,
    table: &impl WideFieldTable,
) -> Vec<Poly512> {
    let xy = x.mul(y, table);
    x.add(y)
        .square(table)
        .mul(&z.square(table), table)
        .add(&xy.mul(z, table))
        .add(&xy.square(table))
        .add(&Sym512::constant(b, table.degree() as usize))
        .coords
}

#[derive(Clone)]
pub struct System512 {
    pub equations: Vec<Poly512>,
    pub n_vars: usize,
    pub summand_bits: usize,
}

/// A one-round degree-four row certificate. Every listed witness occurs in
/// exactly one prolonged row, so that row cannot participate in a
/// combination whose degree-four part vanishes.
#[derive(Clone, Debug)]
pub struct PrivateDegreeFour {
    pub prolonged_rows: usize,
    pub rows_with_degree_four: usize,
    pub degree_four_occurrences: u64,
    pub candidate_columns: usize,
    pub private_witnesses: Vec<(usize, Mono512)>,
    pub unresolved_rows: Vec<usize>,
}

/// A degree-three certificate for rows left by the degree-four pass.
/// Original equations participate in occurrence counts and are retained.
#[derive(Clone, Debug)]
pub struct PrivateDegreeThree {
    pub input_rows: usize,
    pub rows_with_degree_three: usize,
    pub degree_three_occurrences: u64,
    pub candidate_columns: usize,
    pub private_witnesses: Vec<(usize, Mono512)>,
    pub unresolved_rows: Vec<usize>,
}

impl System512 {
    pub fn build(
        basis: &[F2mElement],
        target_x: &F2mElement,
        b: &F2mElement,
        m: usize,
        table: &impl WideFieldTable,
    ) -> Option<Self> {
        let n = table.degree() as usize;
        if !(2..=6).contains(&m)
            || n == 0
            || n > table.max_degree() as usize
            || m.checked_mul(basis.len())?
                .checked_add((m - 2).checked_mul(n)?)?
                > MAX_VARS
        {
            return None;
        }
        let ell = basis.len();
        let xs: Vec<_> = (0..m)
            .map(|i| Sym512::subspace(basis, i * ell, n))
            .collect();
        let inter: Vec<_> = (0..m - 2)
            .map(|i| Sym512::free(m * ell + i * n, n))
            .collect();
        let target = Sym512::constant(target_x, n);
        let mut equations = Vec::with_capacity((m - 1) * n);
        if m == 2 {
            equations.extend(s3(&xs[0], &xs[1], &target, b, table));
        } else {
            equations.extend(s3(&xs[0], &xs[1], &inter[0], b, table));
            for i in 0..m - 3 {
                equations.extend(s3(&inter[i], &xs[i + 2], &inter[i + 1], b, table));
            }
            equations.extend(s3(&inter[m - 3], &xs[m - 1], &target, b, table));
        }
        equations.retain(|p| !p.terms.is_empty());
        Some(Self {
            equations,
            n_vars: m * ell + (m - 2) * n,
            summand_bits: m * ell,
        })
    }

    pub fn all_vanish(&self, assignment: &Mono512) -> bool {
        self.equations.iter().all(|p| !p.eval(assignment))
    }

    pub fn monomial_count(&self) -> usize {
        self.equations.iter().map(|p| p.terms.len()).sum()
    }

    pub fn max_degree(&self) -> u32 {
        self.equations
            .iter()
            .map(Poly512::degree)
            .max()
            .unwrap_or(0)
    }

    /// Fix one complete source-coordinate code while retaining the same
    /// variable layout for the remaining system and its certificates.
    pub fn assign_summand_code(&self, summand: usize, ell: usize, code: u64) -> Option<Self> {
        let start = summand.checked_mul(ell)?;
        if ell == 0
            || ell > 64
            || start.checked_add(ell)? > self.summand_bits
            || (ell < 64 && code >> ell != 0)
        {
            return None;
        }
        let mut zeros = Mono512::default();
        let mut ones = Mono512::default();
        for i in 0..ell {
            let destination = if code >> i & 1 == 1 {
                &mut ones
            } else {
                &mut zeros
            };
            destination.0[(start + i) / 64] |= 1u64 << ((start + i) % 64);
        }
        let mut equations: Vec<_> = self
            .equations
            .iter()
            .map(|p| p.assign_constants(zeros, ones))
            .collect();
        equations.retain(|p| !p.terms.is_empty());
        Some(Self {
            equations,
            n_vars: self.n_vars,
            summand_bits: self.summand_bits,
        })
    }

    /// One selective Macaulay prolongation by original source variables.
    /// Products are formed only from the original equations; newly added
    /// rows are never multiplied again. Boolean idempotence and duplicate
    /// cancellation are handled by `Poly512::mul`.
    pub fn prolongate_source_variables(&self, variables: &[usize]) -> Option<Self> {
        if variables.iter().any(|&v| v >= self.summand_bits) {
            return None;
        }
        let extra = self.equations.len().checked_mul(variables.len())?;
        let mut equations = Vec::with_capacity(self.equations.len().checked_add(extra)?);
        equations.extend_from_slice(&self.equations);
        for &variable in variables {
            let multiplier = Poly512::var(variable);
            for equation in &self.equations {
                let product = equation.mul(&multiplier);
                if !product.terms.is_empty() {
                    equations.push(product);
                }
            }
        }
        Some(Self {
            equations,
            n_vars: self.n_vars,
            summand_bits: self.summand_bits,
        })
    }

    /// Two exact passes over one-round source products. Only a small,
    /// deterministic set of degree-four monomials per row is indexed;
    /// pass two counts each candidate across *all* rows. Missing a
    /// private monomial leaves a row unresolved, never falsely certified.
    pub fn private_degree_four(
        &self,
        variables: &[usize],
        candidates_per_row: usize,
    ) -> Option<PrivateDegreeFour> {
        if candidates_per_row == 0 || variables.iter().any(|&v| v >= self.summand_bits) {
            return None;
        }
        let mut seen_variables = vec![false; self.summand_bits];
        for &variable in variables {
            if std::mem::replace(&mut seen_variables[variable], true) {
                return None;
            }
        }
        let prolonged_rows = self.equations.len().checked_mul(variables.len())?;
        let mut row_candidates = Vec::with_capacity(prolonged_rows);
        let mut candidate_counts: FxMap<Mono512, u32> = FxMap::default();
        let mut rows_with_degree_four = 0;
        let mut degree_four_occurrences = 0u64;
        for &variable in variables {
            let multiplier = Poly512::var(variable);
            for equation in &self.equations {
                let product = equation.mul(&multiplier);
                let degree_four: Vec<_> = product
                    .terms
                    .iter()
                    .copied()
                    .filter(|monomial| monomial.degree() == 4)
                    .collect();
                degree_four_occurrences =
                    degree_four_occurrences.checked_add(degree_four.len() as u64)?;
                rows_with_degree_four += usize::from(!degree_four.is_empty());
                let count = degree_four.len().min(candidates_per_row);
                let mut choices = Vec::with_capacity(count);
                for j in 0..count {
                    let monomial = degree_four[j * degree_four.len() / count];
                    candidate_counts.entry(monomial).or_insert(0);
                    choices.push(monomial);
                }
                row_candidates.push(choices);
            }
        }
        for &variable in variables {
            let multiplier = Poly512::var(variable);
            for equation in &self.equations {
                let product = equation.mul(&multiplier);
                for monomial in product.terms {
                    if let Some(count) = candidate_counts.get_mut(&monomial) {
                        *count = count.checked_add(1)?;
                    }
                }
            }
        }
        let candidate_columns = candidate_counts.len();
        let mut private_witnesses = Vec::new();
        let mut unresolved_rows = Vec::new();
        for (row, choices) in row_candidates.into_iter().enumerate() {
            if let Some(monomial) = choices
                .into_iter()
                .find(|monomial| candidate_counts[monomial] == 1)
            {
                private_witnesses.push((row, monomial));
            } else {
                unresolved_rows.push(row);
            }
        }
        Some(PrivateDegreeFour {
            prolonged_rows,
            rows_with_degree_four,
            degree_four_occurrences,
            candidate_columns,
            private_witnesses,
            unresolved_rows,
        })
    }

    /// Count sampled cubic monomials in the quartic-unresolved products
    /// and all original equations. A private cubic certifies that its
    /// product row is unnecessary for any degree-two-or-lower consequence.
    pub fn private_degree_three_after_quartic(
        &self,
        variables: &[usize],
        quartic_unresolved_rows: &[usize],
        candidates_per_row: usize,
    ) -> Option<PrivateDegreeThree> {
        if candidates_per_row == 0 || variables.iter().any(|&v| v >= self.summand_bits) {
            return None;
        }
        let base_rows = self.equations.len();
        let total = base_rows.checked_mul(variables.len())?;
        let mut previous = None;
        for &row in quartic_unresolved_rows {
            if row >= total || previous.is_some_and(|last| last >= row) {
                return None;
            }
            previous = Some(row);
        }
        let mut row_candidates = Vec::with_capacity(quartic_unresolved_rows.len());
        let mut candidate_counts: FxMap<Mono512, u32> = FxMap::default();
        let mut rows_with_degree_three = 0;
        let mut degree_three_occurrences = 0u64;
        for &row in quartic_unresolved_rows {
            let product =
                self.equations[row % base_rows].mul(&Poly512::var(variables[row / base_rows]));
            let degree_three: Vec<_> = product
                .terms
                .iter()
                .copied()
                .filter(|monomial| monomial.degree() == 3)
                .collect();
            degree_three_occurrences =
                degree_three_occurrences.checked_add(degree_three.len() as u64)?;
            rows_with_degree_three += usize::from(!degree_three.is_empty());
            let count = degree_three.len().min(candidates_per_row);
            let mut choices = Vec::with_capacity(count);
            for j in 0..count {
                let monomial = degree_three[j * degree_three.len() / count];
                candidate_counts.entry(monomial).or_insert(0);
                choices.push(monomial);
            }
            row_candidates.push(choices);
        }
        for equation in &self.equations {
            for &monomial in &equation.terms {
                if monomial.degree() == 3 {
                    if let Some(count) = candidate_counts.get_mut(&monomial) {
                        *count = count.checked_add(1)?;
                    }
                }
            }
        }
        for &row in quartic_unresolved_rows {
            let product =
                self.equations[row % base_rows].mul(&Poly512::var(variables[row / base_rows]));
            for monomial in product.terms {
                if monomial.degree() == 3 {
                    if let Some(count) = candidate_counts.get_mut(&monomial) {
                        *count = count.checked_add(1)?;
                    }
                }
            }
        }
        let candidate_columns = candidate_counts.len();
        let mut private_witnesses = Vec::new();
        let mut unresolved_rows = Vec::new();
        for (&row, choices) in quartic_unresolved_rows.iter().zip(row_candidates) {
            if let Some(monomial) = choices
                .into_iter()
                .find(|monomial| candidate_counts[monomial] == 1)
            {
                private_witnesses.push((row, monomial));
            } else {
                unresolved_rows.push(row);
            }
        }
        Some(PrivateDegreeThree {
            input_rows: quartic_unresolved_rows.len(),
            rows_with_degree_three,
            degree_three_occurrences,
            candidate_columns,
            private_witnesses,
            unresolved_rows,
        })
    }

    /// Keep every original equation and only the prolonged rows that lack a
    /// private degree-four certificate. This preserves all consequences of
    /// degree at most three from the complete one-round row span.
    pub fn source_prolongation_core(
        &self,
        variables: &[usize],
        unresolved_rows: &[usize],
    ) -> Option<Self> {
        if variables.iter().any(|&v| v >= self.summand_bits) {
            return None;
        }
        let base_rows = self.equations.len();
        let total = base_rows.checked_mul(variables.len())?;
        let mut equations = self.equations.clone();
        let mut previous = None;
        for &row in unresolved_rows {
            if row >= total || previous.is_some_and(|last| last >= row) {
                return None;
            }
            previous = Some(row);
            let product =
                self.equations[row % base_rows].mul(&Poly512::var(variables[row / base_rows]));
            if !product.terms.is_empty() {
                equations.push(product);
            }
        }
        Some(Self {
            equations,
            n_vars: self.n_vars,
            summand_bits: self.summand_bits,
        })
    }

    /// One own-degree Macaulay reduction; a column-cap stop is inconclusive.
    pub fn root_reduce(&self) -> RootReduction {
        self.root_reduce_with_column_cap(MAX_ROOT_COLS)
    }

    /// The same exact reduction with an explicit resource cap. The cap only
    /// decides whether to attempt the matrix; it does not truncate columns.
    pub fn root_reduce_with_column_cap(&self, column_cap: usize) -> RootReduction {
        let mut index: FxMap<Mono512, usize> = FxMap::default();
        for p in &self.equations {
            for &m in &p.terms {
                index.entry(m).or_insert(0);
                if index.len() > column_cap {
                    return RootReduction::ColumnLimit {
                        columns: index.len(),
                    };
                }
            }
        }
        let mut columns: Vec<_> = index.keys().copied().collect();
        columns.sort_unstable_by(|a, b| b.degree().cmp(&a.degree()).then(a.cmp(b)));
        for (i, &m) in columns.iter().enumerate() {
            index.insert(m, i);
        }
        let words = columns.len().div_ceil(64);
        let mut matrix: Vec<Vec<u64>> = self
            .equations
            .iter()
            .map(|p| {
                let mut row = vec![0; words];
                for m in &p.terms {
                    let c = index[m];
                    row[c / 64] |= 1u64 << (c % 64);
                }
                row
            })
            .collect();
        let mut ops = 0;
        // Echelon form is enough: all columns following a degree-one
        // pivot have degree at most one, and a constant pivot is a
        // contradiction regardless of back-substitution.
        let rank = gf2_elim::eliminate(
            &mut matrix,
            columns.len(),
            false,
            gf2_elim::Config::from_env(),
            &mut ops,
        );
        let mut linear_equations = Vec::new();
        let mut contradiction = false;
        for row in &matrix[..rank] {
            let lead = row
                .iter()
                .enumerate()
                .find(|(_, w)| **w != 0)
                .map(|(w, bits)| w * 64 + bits.trailing_zeros() as usize)
                .expect("nonzero pivot");
            if columns[lead].degree() <= 1 {
                let mut variables = Vec::new();
                let mut constant = false;
                for (w, &word) in row.iter().enumerate() {
                    let mut bits = word;
                    while bits != 0 {
                        let c = w * 64 + bits.trailing_zeros() as usize;
                        bits &= bits - 1;
                        let mono = columns[c];
                        match mono.degree() {
                            0 => constant = true,
                            1 => {
                                let (word_index, &value) = mono
                                    .0
                                    .iter()
                                    .enumerate()
                                    .find(|(_, value)| **value != 0)
                                    .expect("linear monomial");
                                variables.push(word_index * 64 + value.trailing_zeros() as usize);
                            }
                            _ => unreachable!("degree-one pivot cannot have a higher-degree tail"),
                        }
                    }
                }
                variables.sort_unstable();
                contradiction |= variables.is_empty() && constant;
                linear_equations.push(LinearEquation {
                    variables,
                    constant,
                });
            }
        }
        RootReduction::Reduced {
            columns: columns.len(),
            rank,
            linear: linear_equations.len(),
            linear_equations,
            contradiction,
            xor_ops: ops,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearEquation {
    /// XOR of these variable bits equals `constant`.
    pub variables: Vec<usize>,
    pub constant: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RootReduction {
    ColumnLimit {
        columns: usize,
    },
    Reduced {
        columns: usize,
        rank: usize,
        linear: usize,
        linear_equations: Vec<LinearEquation>,
        contradiction: bool,
        xor_ops: u64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn word_boundaries_and_boolean_cancellation() {
        for v in [0, 63, 64, 127, 128, 383, 384, 427, 511] {
            let x = Poly512::var(v);
            assert_eq!(x.mul(&x), x);
            assert_eq!(x.add(&x), Poly512::zero());
            let mut assignment = Mono512::default();
            assignment.0[v / 64] |= 1 << (v % 64);
            assert!(x.eval(&assignment));
            assert!(!x.eval(&Mono512::default()));
            assert_eq!(x.mul(&Poly512::one()), x);
        }
        let x = Poly512::var(64);
        let y = Poly512::var(384);
        let product = x.mul(&y);
        assert_eq!(product.terms[0].degree(), 2);
        assert_eq!(product.add(&product), Poly512::zero());
    }

    #[test]
    fn exact_linear_consequence_replays() {
        let system = System512 {
            equations: vec![Poly512::var(64)
                .add(&Poly512::var(384))
                .add(&Poly512::one())],
            n_vars: 512,
            summand_bits: 96,
        };
        let RootReduction::Reduced {
            rank,
            linear_equations,
            contradiction,
            ..
        } = system.root_reduce()
        else {
            panic!("small matrix must reduce")
        };
        assert_eq!(rank, 1);
        assert!(!contradiction);
        assert_eq!(
            linear_equations,
            vec![LinearEquation {
                variables: vec![64, 384],
                constant: true
            }]
        );
    }

    #[test]
    fn source_code_substitution_preserves_evaluation_and_cancellation() {
        let p = Poly512::var(63)
            .mul(&Poly512::var(64))
            .add(&Poly512::var(64))
            .add(&Poly512::var(384));
        let mut zeros = Mono512::default();
        zeros.0[0] = 1u64 << 63;
        let reduced = p.assign_constants(zeros, Mono512::default());
        assert_eq!(reduced, Poly512::var(64).add(&Poly512::var(384)));
        let mut ones = Mono512::default();
        ones.0[0] = 1u64 << 63;
        assert_eq!(
            p.assign_constants(Mono512::default(), ones),
            Poly512::var(384)
        );
        let system = System512 {
            equations: vec![p],
            n_vars: 512,
            summand_bits: 96,
        };
        let fixed = system
            .assign_summand_code(3, 16, 0)
            .expect("v63 is in summand 3");
        assert_eq!(fixed.equations[0], reduced);
    }

    #[test]
    fn selective_source_prolongation_preserves_planted_zero() {
        let original = System512 {
            equations: vec![Poly512::var(0).add(&Poly512::var(1)).add(&Poly512::one())],
            n_vars: 3,
            summand_bits: 2,
        };
        let mut assignment = Mono512::default();
        assignment.0[0] = 1;
        assert!(original.all_vanish(&assignment));
        let prolonged = original.prolongate_source_variables(&[0, 1]).unwrap();
        assert_eq!(prolonged.equations.len(), 3);
        assert!(prolonged.all_vanish(&assignment));
        assert!(original.prolongate_source_variables(&[2]).is_none());
    }

    #[test]
    fn private_degree_four_distinguishes_unique_and_shared_columns() {
        let cubic = Poly512::var(0).mul(&Poly512::var(1)).mul(&Poly512::var(2));
        let unique = System512 {
            equations: vec![cubic.clone()],
            n_vars: 5,
            summand_bits: 5,
        };
        let certificate = unique.private_degree_four(&[3], 32).unwrap();
        assert_eq!(certificate.prolonged_rows, 1);
        assert_eq!(certificate.private_witnesses.len(), 1);
        assert!(certificate.unresolved_rows.is_empty());
        let shared = System512 {
            equations: vec![cubic.clone(), cubic],
            n_vars: 5,
            summand_bits: 5,
        };
        let certificate = shared.private_degree_four(&[3], 32).unwrap();
        assert!(certificate.private_witnesses.is_empty());
        assert_eq!(certificate.unresolved_rows, vec![0, 1]);
        let core = shared
            .source_prolongation_core(&[3], &certificate.unresolved_rows)
            .unwrap();
        assert_eq!(core.equations.len(), 4);
        assert!(shared.source_prolongation_core(&[3], &[1, 0]).is_none());
    }

    #[test]
    fn private_degree_three_counts_original_equations() {
        let quadratic = Poly512::var(0).mul(&Poly512::var(1));
        let cubic = quadratic.mul(&Poly512::var(2));
        let unique = System512 {
            equations: vec![quadratic.clone()],
            n_vars: 3,
            summand_bits: 3,
        };
        let certificate = unique
            .private_degree_three_after_quartic(&[2], &[0], 32)
            .unwrap();
        assert_eq!(certificate.private_witnesses, vec![(0, cubic.terms[0])]);
        assert!(certificate.unresolved_rows.is_empty());

        let shared = System512 {
            equations: vec![quadratic, cubic],
            n_vars: 3,
            summand_bits: 3,
        };
        let certificate = shared
            .private_degree_three_after_quartic(&[2], &[0, 1], 32)
            .unwrap();
        assert!(certificate.private_witnesses.is_empty());
        assert_eq!(certificate.unresolved_rows, vec![0, 1]);
    }
}
