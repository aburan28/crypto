// Diagnostic access only: appended to an immutable historical source copy.
// The audited production implementations above are unchanged.
pub fn diagnostic_production_decisive(
    polys: &[F2BoolPoly],
    n_vars: usize,
    degree: u32,
    criterion: RowCriterion,
) -> Option<Vec<F2BoolPoly>> {
    matrix_f4_f2_solver_consequences(polys, n_vars, degree, criterion)
}

pub fn diagnostic_production_substitute(p: &F2BoolPoly, variable: u32, value: bool) -> F2BoolPoly {
    substitute_in_solve(p, variable, value, true)
}
