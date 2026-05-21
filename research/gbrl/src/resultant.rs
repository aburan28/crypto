// Sylvester resultant. Eliminates a single variable from two polynomials
// living in the same ambient (nvars unchanged in the result; the eliminated
// variable has exponent 0 in every output term).
//
// Used to chain Semaev S_m polynomials: S_{m+1} = Res_Y(S_{m-k}(..., Y),
// S_{k+2}(..., Y)). Cofactor expansion on the Sylvester matrix; cubic in
// matrix size and exponential in nesting, so fine for small m (4, 5, maybe 6)
// and rapidly intractable beyond that — at which point you want subresultant
// PRS or modular interpolation.

use crate::field::Fp;
use crate::monomial::Monomial;
use crate::poly::{Poly, Term};

// View `p` as a univariate polynomial in variable `y_idx`. Returns
// `coefs[d]` = the coefficient polynomial of `y_idx^d`, in the same ambient
// (nvars unchanged) but with exponent 0 on `y_idx`.
fn coefs_in_var(p: &Poly, y_idx: usize) -> Vec<Poly> {
    let max_deg = p.terms.iter().map(|t| t.mono.exp[y_idx]).max().unwrap_or(0) as usize;
    let mut bucketed: Vec<Vec<Term>> = vec![Vec::new(); max_deg + 1];
    for t in &p.terms {
        let d = t.mono.exp[y_idx] as usize;
        let mut new_exp = t.mono.exp.clone();
        new_exp[y_idx] = 0;
        bucketed[d].push(Term { coef: t.coef, mono: Monomial { exp: new_exp } });
    }
    bucketed.into_iter().map(|ts| Poly::from_terms(ts, p.nvars)).collect()
}

// Build the (m+n) × (m+n) Sylvester matrix from coefficient lists.
// `p_coefs[i]` is the coef of Y^i in p, similarly for q.
fn sylvester(p_coefs: &[Poly], q_coefs: &[Poly]) -> Vec<Vec<Poly>> {
    let m = p_coefs.len().saturating_sub(1);  // deg_Y p
    let n = q_coefs.len().saturating_sub(1);  // deg_Y q
    let size = m + n;
    let nvars = p_coefs[0].nvars;
    let mut mat = vec![vec![Poly::zero(nvars); size]; size];
    // Top n rows: shifted copies of p's coefficients (highest-degree first).
    for i in 0..n {
        for j in 0..=m {
            mat[i][i + j] = p_coefs[m - j].clone();
        }
    }
    // Bottom m rows: shifted copies of q's coefficients.
    for i in 0..m {
        for j in 0..=n {
            mat[n + i][i + j] = q_coefs[n - j].clone();
        }
    }
    mat
}

// Determinant by cofactor expansion along row 0. O(n!) which is fine for
// n ≤ ~5; switch to Bareiss division-free reduction for bigger matrices.
fn det(mat: &[Vec<Poly>]) -> Poly {
    let n = mat.len();
    if n == 0 {
        // det of empty matrix = 1; placeholder ambient (1-var) is fine here.
        return Poly::from_terms(
            vec![Term { coef: Fp::one(), mono: Monomial::one(1) }],
            1,
        );
    }
    if n == 1 {
        return mat[0][0].clone();
    }
    let nvars = mat[0][0].nvars;
    let mut acc = Poly::zero(nvars);
    for j in 0..n {
        if mat[0][j].is_zero() { continue; }
        let minor: Vec<Vec<Poly>> = mat[1..].iter().map(|row| {
            row.iter().enumerate()
                .filter(|(c, _)| *c != j)
                .map(|(_, p)| p.clone())
                .collect()
        }).collect();
        let sub = det(&minor);
        let term = mat[0][j].mul(&sub);
        if j % 2 == 0 { acc = acc.add(&term); } else { acc = acc.sub(&term); }
    }
    acc
}

// Sylvester resultant Res_{y_idx}(p, q). Output has nvars = p.nvars (same
// ambient as inputs); the eliminated variable always has exponent 0.
pub fn resultant(p: &Poly, q: &Poly, y_idx: usize) -> Poly {
    assert_eq!(p.nvars, q.nvars);
    let pc = coefs_in_var(p, y_idx);
    let qc = coefs_in_var(q, y_idx);
    let mat = sylvester(&pc, &qc);
    det(&mat)
}
