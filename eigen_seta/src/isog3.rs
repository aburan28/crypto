//! Chains of 3-isogenies between Montgomery curves, x-only, constant time
//! (fixed public strategy, branch-free formulas).

use crate::field::Field;
use crate::mont::{eval_3_isog, get_3_isog, xtpl_iter, ProjCurve, XPoint};
use crate::strategy::{optimal, traverse};

/// Codomain of the 3^n-isogeny with kernel <K> (K of order 3^n), and the
/// images of `pts`.
pub fn chain3<F: Field>(curve: &ProjCurve<F>, k: &XPoint<F>, n: usize, pts: &mut [XPoint<F>]) -> ProjCurve<F> {
    let strategy = optimal(n, 1.6, 1.0); // tripling ~ 1.6 evaluations
    let mut cur = *curve;
    let curs = core::cell::RefCell::new(Vec::<ProjCurve<F>>::new());
    // The multiplication closure needs the curve of the current level: track
    // it by recording codomains as isogenies are built.
    let _ = traverse(
        n,
        &strategy,
        *k,
        pts,
        |x, s| {
            let c = curs.borrow().last().copied().unwrap_or(cur);
            xtpl_iter(x, s, &c)
        },
        |x| {
            let (cod, coeff) = get_3_isog(x);
            curs.borrow_mut().push(cod);
            coeff
        },
        |coeff, x| eval_3_isog(x, coeff),
    );
    if let Some(c) = curs.borrow().last() {
        cur = *c;
    }
    cur
}
