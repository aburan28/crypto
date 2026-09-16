//! First-fall degree / Macaulay rank profiles over `F_p`.
//!
//! Operational definition matches the Boolean harness in
//! `koblitz_groebner::first_fall_degree`: smallest `D ≥ 2` whose Macaulay
//! matrix has `rank < rows` and `rank < cols` (non-trivial syzygy before
//! column saturation).

use crate::field::Fp;
use crate::monomial::Monomial;
use crate::poly::Poly;
use std::collections::HashMap;

const MAX_MACAULAY_ROWS: usize = 8_000;
const MAX_MACAULAY_COLS: usize = 4_000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MacaulayProfile {
    pub degree: u32,
    pub rows: usize,
    pub cols: usize,
    pub rank: usize,
}

impl MacaulayProfile {
    pub fn syzygies(&self) -> usize {
        self.rows.saturating_sub(self.rank)
    }
}

/// Enumerate all monomials of total degree ≤ `max_deg` in `nvars` variables
/// (including the constant 1).
pub fn monomials_up_to(nvars: usize, max_deg: u32) -> Vec<Monomial> {
    let mut out = Vec::new();
    let mut exp = vec![0u32; nvars];
    fn rec(out: &mut Vec<Monomial>, exp: &mut [u32], i: usize, remaining: u32) {
        if i == exp.len() {
            out.push(Monomial { exp: exp.to_vec() });
            return;
        }
        for e in 0..=remaining {
            exp[i] = e;
            rec(out, exp, i + 1, remaining - e);
        }
        exp[i] = 0;
    }
    rec(&mut out, &mut exp, 0, max_deg);
    out
}

fn build_macaulay(
    polys: &[Poly],
    nvars: usize,
    degree: u32,
) -> Option<(Vec<Monomial>, Vec<Vec<Fp>>)> {
    let mut row_polys: Vec<Poly> = Vec::new();
    for p in polys {
        let pdeg = p.total_degree();
        if pdeg > degree {
            continue;
        }
        for mult in monomials_up_to(nvars, degree - pdeg) {
            let shifted = p.mul_term(Fp::one(), &mult);
            if !shifted.is_zero() {
                row_polys.push(shifted);
            }
            if row_polys.len() > MAX_MACAULAY_ROWS {
                return None;
            }
        }
    }
    if row_polys.is_empty() {
        return Some((Vec::new(), Vec::new()));
    }

    let mut cols: Vec<Monomial> = row_polys
        .iter()
        .flat_map(|p| p.terms.iter().map(|t| t.mono.clone()))
        .collect();
    cols.sort_by(|a, b| b.cmp(a));
    cols.dedup();
    if cols.len() > MAX_MACAULAY_COLS {
        return None;
    }

    let index: HashMap<Monomial, usize> = cols
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, m)| (m, i))
        .collect();
    let matrix: Vec<Vec<Fp>> = row_polys
        .iter()
        .map(|p| {
            let mut row = vec![Fp::zero(); cols.len()];
            for term in &p.terms {
                row[index[&term.mono]] = term.coef;
            }
            row
        })
        .collect();
    Some((cols, matrix))
}

fn rank_fp(matrix: &mut [Vec<Fp>], ncols: usize) -> usize {
    let mut pivot_row = 0usize;
    for col in 0..ncols {
        let Some(found) = (pivot_row..matrix.len()).find(|&r| !matrix[r][col].is_zero()) else {
            continue;
        };
        matrix.swap(pivot_row, found);
        let inv = matrix[pivot_row][col].inv().expect("nonzero pivot");
        for c in col..ncols {
            matrix[pivot_row][c] = matrix[pivot_row][c] * inv;
        }
        for r in 0..matrix.len() {
            if r == pivot_row || matrix[r][col].is_zero() {
                continue;
            }
            let factor = matrix[r][col];
            for c in col..ncols {
                matrix[r][c] = matrix[r][c] - factor * matrix[pivot_row][c];
            }
        }
        pivot_row += 1;
        if pivot_row == matrix.len() {
            break;
        }
    }
    pivot_row
}

pub fn macaulay_profile(polys: &[Poly], nvars: usize, degree: u32) -> Option<MacaulayProfile> {
    let (cols, mut matrix) = build_macaulay(polys, nvars, degree)?;
    let rows = matrix.len();
    let rank = if matrix.is_empty() {
        0
    } else {
        rank_fp(&mut matrix, cols.len())
    };
    Some(MacaulayProfile {
        degree,
        rows,
        cols: cols.len(),
        rank,
    })
}

/// First fall degree of `polys` over `F_p`, scanning `D = 2..=d_max`.
pub fn first_fall_degree(
    polys: &[Poly],
    nvars: usize,
    d_max: u32,
) -> (Option<u32>, Vec<MacaulayProfile>) {
    let mut fall = None;
    let mut profiles = Vec::new();
    for d in 2..=d_max {
        let Some(prof) = macaulay_profile(polys, nvars, d) else {
            break;
        };
        if fall.is_none() && prof.rank < prof.rows && prof.rank < prof.cols {
            fall = Some(d);
        }
        profiles.push(prof);
    }
    (fall, profiles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Term;

    #[test]
    fn monomials_up_to_two_vars_deg1() {
        let ms = monomials_up_to(2, 1);
        assert_eq!(ms.len(), 3); // 1, x, y
    }

    #[test]
    fn ffd_detects_duplicate_generator_syzygy() {
        // Two identical linear generators → row dependency (syzygies > 0).
        // FFD may be None when occurring columns saturate (rank == cols).
        let x = Monomial { exp: vec![1, 0] };
        let polys = vec![
            Poly::from_terms(vec![Term { coef: Fp::one(), mono: x.clone() }], 2),
            Poly::from_terms(vec![Term { coef: Fp::one(), mono: x }], 2),
        ];
        let (_fall, profiles) = first_fall_degree(&polys, 2, 4);
        assert!(!profiles.is_empty());
        assert!(profiles.iter().any(|p| p.syzygies() > 0));
    }
}
