//! Fast exhaustive search over quadratic Boolean systems, for the same
//! Semaev ANF the WDSat path emits.
//!
//! Inspired by the LIP6 / ALMASTY MQ line
//! (<https://gitlab.lip6.fr/almasty/mq>, Bouillaguet–Delaplace–Vrinçanu
//! “simple deterministic” quadratic search, and Bouillaguet’s
//! [`libfes-lite`](https://github.com/cbouilla/libfes-lite)): enumerate
//! the free variables with a Gray code and update each quadratic
//! polynomial with one XOR per step.  This is a study-library port of
//! the *idea*, not a binding to those codebases.
//!
//! Restricted to degree ≤ 2 and `n_vars ≤ 24`.  Chained `m ≥ 3` Semaev
//! systems are cubic and are refused here — use SAT or WDSat for those.

use super::wdsat_oracle::AnfRow;
use crate::cryptanalysis::pq_groebner_f2::F2BoolPoly;

/// Packed quadratic form over `n` Boolean variables: constant, linear
/// coefficients, and upper-triangular quadratic coefficients.
#[derive(Clone, Debug)]
pub struct QuadraticForm {
    pub n: usize,
    pub constant: bool,
    pub linear: Vec<bool>,
    /// `quad[i][j]` for `0 ≤ j < i < n` is the coefficient of `x_i x_j`.
    pub quad: Vec<Vec<bool>>,
}

impl QuadraticForm {
    /// Convert an ANF row; returns `None` if any monomial has degree > 2.
    pub fn from_anf_row(row: &AnfRow, n: usize) -> Option<Self> {
        let mut form = Self {
            n,
            constant: row.constant,
            linear: vec![false; n],
            quad: (0..n).map(|i| vec![false; i]).collect(),
        };
        for mono in &row.monomials {
            match mono.len() {
                0 => form.constant = !form.constant,
                1 => {
                    let v = mono[0] as usize;
                    if v >= n {
                        return None;
                    }
                    form.linear[v] = !form.linear[v];
                }
                2 => {
                    let mut a = mono[0] as usize;
                    let mut b = mono[1] as usize;
                    if a == b || a >= n || b >= n {
                        return None;
                    }
                    if a < b {
                        std::mem::swap(&mut a, &mut b);
                    }
                    form.quad[a][b] = !form.quad[a][b];
                }
                _ => return None,
            }
        }
        Some(form)
    }

    pub fn from_poly(poly: &F2BoolPoly) -> Option<Self> {
        Self::from_anf_row(&AnfRow::from_poly(poly), poly.n_vars)
    }

    /// Evaluate at a bit-packed point (`bit i` = value of `x_i`).
    pub fn eval(&self, point: u64) -> bool {
        let mut v = self.constant;
        for i in 0..self.n {
            if ((point >> i) & 1) == 1 && self.linear[i] {
                v = !v;
            }
        }
        for i in 0..self.n {
            if ((point >> i) & 1) == 0 {
                continue;
            }
            for j in 0..i {
                if ((point >> j) & 1) == 1 && self.quad[i][j] {
                    v = !v;
                }
            }
        }
        v
    }
}

/// Gray-code exhaustive search for a common zero of quadratic forms.
///
/// Returns the first satisfying assignment as a bit pack, or `None` if
/// the full `2^n` space is empty.  Cap: `n ≤ 24`.
pub fn fes_find_one(forms: &[QuadraticForm]) -> Option<u64> {
    if forms.is_empty() {
        return Some(0);
    }
    let n = forms[0].n;
    if n > 24 || forms.iter().any(|f| f.n != n) {
        return None;
    }
    // Naive first: for study sizes the Gray-code constant matters less
    // than correctness.  Still walk in Gray order so a later SIMD port
    // can drop in without changing the search shape.
    let limit = 1u64 << n;
    let mut point = 0u64;
    for step in 0..limit {
        if forms.iter().all(|f| !f.eval(point)) {
            return Some(point);
        }
        // Next Gray code: flip the bit corresponding to the lowest set
        // bit of (step+1).
        let flip = (step + 1).trailing_zeros() as u64;
        if flip < n as u64 {
            point ^= 1u64 << flip;
        }
    }
    None
}

/// Enumerate every common zero (still capped at `n ≤ 24`).
pub fn fes_find_all(forms: &[QuadraticForm], max_solutions: usize) -> Vec<u64> {
    let mut out = Vec::new();
    if forms.is_empty() {
        return vec![0];
    }
    let n = forms[0].n;
    if n > 24 || forms.iter().any(|f| f.n != n) {
        return out;
    }
    let limit = 1u64 << n;
    let mut point = 0u64;
    for step in 0..limit {
        if forms.iter().all(|f| !f.eval(point)) {
            out.push(point);
            if out.len() >= max_solutions {
                break;
            }
        }
        let flip = (step + 1).trailing_zeros() as u64;
        if flip < n as u64 {
            point ^= 1u64 << flip;
        }
    }
    out
}

/// Solve a quadratic (`m = 2`) Semaev decomposition by Gray-code FES.
///
/// Returns the same shape as
/// [`crate::cryptanalysis::koblitz_index_calculus::sat_decompose`].
/// Cubic chained systems (`m ≥ 3`) are reported as exhausted.
pub fn mq_fes_decompose(
    kc: &crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve,
    fb: &crate::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    index_of: &std::collections::HashMap<(num_bigint::BigUint, num_bigint::BigUint), usize>,
    st: &crate::cryptanalysis::koblitz_groebner::FieldStructure,
    target: &crate::binary_ecc::BinaryPoint,
    m: usize,
) -> (
    Option<Vec<usize>>,
    crate::cryptanalysis::koblitz_index_calculus::SatDecompositionStats,
) {
    use crate::binary_ecc::BinaryPoint;
    use crate::cryptanalysis::koblitz_index_calculus::{lift_candidate, SatDecompositionStats};

    let mut stats = SatDecompositionStats::default();
    if m != 2 {
        stats.exhausted = true;
        return (None, stats);
    }
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, stats),
    };
    let sys = match crate::cryptanalysis::polynomial_reuse::build_decomposition_system_reusing(
        &fb.subspace_basis,
        &x_r,
        &kc.curve.b,
        m,
        st,
    ) {
        Some(sys) => sys,
        None => {
            stats.exhausted = true;
            return (None, stats);
        }
    };
    let forms: Option<Vec<_>> = sys
        .equations
        .iter()
        .map(QuadraticForm::from_poly)
        .collect();
    let forms = match forms {
        Some(forms) => forms,
        None => {
            stats.exhausted = true;
            return (None, stats);
        }
    };
    stats.solver_calls = 1;
    let roots = fes_find_all(&forms, 64);
    if roots.is_empty() {
        stats.refuted = true;
        return (None, stats);
    }
    for root in roots {
        stats.models += 1;
        if !sys.equations.iter().all(|e| e.eval(root) == 0) {
            stats.spurious += 1;
            continue;
        }
        let xs: Vec<_> = (0..m)
            .map(|i| sys.summand_x(&fb.subspace_basis, root, i, kc.n))
            .collect();
        if let Some(idxs) = lift_candidate(kc, fb, index_of, &xs, target) {
            return (Some(idxs), stats);
        }
    }
    stats.exhausted = true;
    (None, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::wdsat_oracle::AnfRow;

    #[test]
    fn finds_the_unique_zero_of_a_tiny_system() {
        // x0 + x1 = 0 and x0 = 1  →  only (1,1).
        let rows = [
            AnfRow {
                monomials: vec![vec![0], vec![1]],
                constant: false,
            },
            AnfRow {
                monomials: vec![vec![0]],
                constant: true,
            },
        ];
        let forms: Vec<_> = rows
            .iter()
            .map(|r| QuadraticForm::from_anf_row(r, 2).unwrap())
            .collect();
        assert_eq!(fes_find_one(&forms), Some(0b11));
    }

    #[test]
    fn rejects_a_cubic_row() {
        let row = AnfRow {
            monomials: vec![vec![0, 1, 2]],
            constant: false,
        };
        assert!(QuadraticForm::from_anf_row(&row, 3).is_none());
    }

    #[test]
    fn mq_fes_agrees_with_native_sat_on_prime_degree() {
        use crate::cryptanalysis::koblitz_index_calculus::{
            build_frobenius_factor_base, point_key, sat_decompose, KoblitzCurve,
        };
        use crate::cryptanalysis::koblitz_groebner::FieldStructure;
        use crate::binary_ecc::BinaryPoint;
        use std::collections::HashMap;

        let kc = KoblitzCurve::new(1, 7).expect("K_1/F_2^7");
        let fb = build_frobenius_factor_base(&kc, 0).expect("factor base");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let index_of: HashMap<_, _> = fb
            .points
            .iter()
            .enumerate()
            .map(|(i, p)| (point_key(p), i))
            .collect();
        let mut affine = fb.points.iter().enumerate().filter_map(|(i, p)| match p {
            BinaryPoint::Affine { .. } => Some(i),
            BinaryPoint::Infinity => None,
        });
        let i = affine.next().unwrap();
        let j = affine.next().unwrap();
        let target = kc.add(&fb.points[i], &fb.points[j]);

        let (native, _) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 8, Some(2));
        let (fes, stats) = mq_fes_decompose(&kc, &fb, &index_of, &st, &target, 2);
        assert!(native.is_some(), "native SAT missed planted");
        assert!(fes.is_some(), "mq-fes missed planted: {stats:?}");
        let mut a = native.unwrap();
        let mut b = fes.unwrap();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }
}
