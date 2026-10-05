//! Fixed-width Boolean system for the direct n83 six-summand S3 chain.
//! This admission path is independent of the 128-variable wide solver.

use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::fx_hash::FxMap;
use crate::cryptanalysis::gf2_elim;
use crate::cryptanalysis::wide_groebner::WideFieldTable;

pub const MAX_VARS: usize = 512;
pub const MAX_ROOT_COLS: usize = 300_000;

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

    /// One own-degree Macaulay reduction; a column-cap stop is inconclusive.
    pub fn root_reduce(&self) -> RootReduction {
        let mut index: FxMap<Mono512, usize> = FxMap::default();
        for p in &self.equations {
            for &m in &p.terms {
                index.entry(m).or_insert(0);
                if index.len() > MAX_ROOT_COLS {
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
        let mut linear = 0;
        let mut contradiction = false;
        for row in &matrix[..rank] {
            let lead = row
                .iter()
                .enumerate()
                .find(|(_, w)| **w != 0)
                .map(|(w, bits)| w * 64 + bits.trailing_zeros() as usize)
                .expect("nonzero pivot");
            if columns[lead].degree() <= 1 {
                linear += 1;
                if columns[lead].degree() == 0 {
                    contradiction = true;
                }
            }
        }
        RootReduction::Reduced {
            columns: columns.len(),
            rank,
            linear,
            contradiction,
            xor_ops: ops,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RootReduction {
    ColumnLimit {
        columns: usize,
    },
    Reduced {
        columns: usize,
        rank: usize,
        linear: usize,
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
}
