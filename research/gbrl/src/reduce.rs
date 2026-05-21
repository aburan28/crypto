use crate::poly::{Poly, Term};

// Multivariate division. Top-reduce f modulo the basis, returning (remainder, op_count).
// Op count = number of leading-term cancellations performed (a coarse proxy for arith cost).
pub fn reduce(f: &Poly, basis: &[Poly]) -> (Poly, usize) {
    let mut p = f.clone();
    let mut r_terms: Vec<Term> = Vec::new();
    let mut steps = 0;

    loop {
        if p.is_zero() { break; }
        let lt = p.lt().unwrap().clone();

        let mut reduced = false;
        for g in basis {
            if g.is_zero() { continue; }
            let gt = g.lt().unwrap();
            if let Some(m) = lt.mono.div(&gt.mono) {
                // p -= (lt.coef / gt.coef) * m * g
                let c = lt.coef * gt.coef.inv().expect("nonzero leading coef");
                let sub = g.mul_term(c, &m);
                p = p.sub(&sub);
                steps += 1;
                reduced = true;
                break;
            }
        }

        if !reduced {
            // LT(p) cannot be reduced — move it to the remainder. Because monomials
            // strictly decrease across iterations, r_terms stays sorted descending.
            r_terms.push(lt);
            p.terms.remove(0);
        }
    }

    (Poly { terms: r_terms, nvars: f.nvars }, steps)
}
