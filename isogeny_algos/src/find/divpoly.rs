//! Kernel finding by factoring the division polynomial psi_l (Schoof / Elkies / Couveignes
//! l-torsion style): every F-rational cyclic subgroup of order l has kernel polynomial
//! equal to a product of irreducible factors of psi_l of total degree (l-1)/2.
use crate::curve::*;
use crate::field::{Field, Rng};
use crate::kernel::kohel::kohel;
use crate::poly::{self, Poly};
use std::collections::HashMap;

/// g_n(x): psi_n for odd n, psi_n/(2y) for even n, for y^2 = x^3+ax+b.
pub fn division_poly<F: Field>(f: &F, e: &Curve<F::E>, ell: usize) -> Poly<F> {
    let c = |v: u64| f.from_u64(v);
    let fx: Poly<F> = vec![e.b, e.a, f.zero(), f.one()];
    let f16: Poly<F> = poly::scale(f, &poly::mul(f, &fx, &fx), c(16));
    let a2 = f.mul(e.a, e.a);
    let mut memo: HashMap<usize, Poly<F>> = HashMap::new();
    memo.insert(0, vec![]);
    memo.insert(1, vec![f.one()]);
    memo.insert(2, vec![f.one()]);
    memo.insert(
        3,
        vec![
            f.neg(a2),
            f.mul(c(12), e.b),
            f.mul(c(6), e.a),
            f.zero(),
            c(3),
        ],
    );
    let a3 = f.mul(a2, e.a);
    let b2 = f.mul(e.b, e.b);
    memo.insert(
        4,
        poly::scale(
            f,
            &vec![
                f.sub(f.neg(f.mul(c(8), b2)), a3),
                f.neg(f.mul(c(4), f.mul(e.a, e.b))),
                f.neg(f.mul(c(5), a2)),
                f.mul(c(20), e.b),
                f.mul(c(5), e.a),
                f.zero(),
                f.one(),
            ],
            c(2),
        ),
    );
    // psi_4 = 4y(...) => g_4 = 2*(...)
    fn get<F: Field>(
        f: &F,
        n: usize,
        memo: &mut HashMap<usize, Poly<F>>,
        f16: &Poly<F>,
    ) -> Poly<F> {
        if let Some(v) = memo.get(&n) {
            return v.clone();
        }
        let r = if n % 2 == 1 {
            let m = (n - 1) / 2;
            let (a, b) = (get(f, m + 2, memo, f16), get(f, m, memo, f16));
            let (c1, d1) = (get(f, m - 1, memo, f16), get(f, m + 1, memo, f16));
            let cube = |p: &Poly<F>| poly::mul(f, p, &poly::mul(f, p, p));
            let t1 = poly::mul(f, &a, &cube(&b));
            let t2 = poly::mul(f, &c1, &cube(&d1));
            if m.is_multiple_of(2) {
                poly::sub(f, &poly::mul(f, f16, &t1), &t2)
            } else {
                poly::sub(f, &t1, &poly::mul(f, f16, &t2))
            }
        } else {
            let m = n / 2;
            let gm = get(f, m, memo, f16);
            let (gp2, gm1) = (get(f, m + 2, memo, f16), get(f, m - 1, memo, f16));
            let (gm2, gp1) = (get(f, m - 2, memo, f16), get(f, m + 1, memo, f16));
            let l = poly::mul(f, &gp2, &poly::mul(f, &gm1, &gm1));
            let r = poly::mul(f, &gm2, &poly::mul(f, &gp1, &gp1));
            poly::mul(f, &gm, &poly::sub(f, &l, &r))
        };
        memo.insert(n, r.clone());
        r
    }
    get(f, ell, &mut memo, &f16)
}

/// All F-rational kernel polynomials of degree-l isogenies (l odd prime).
pub fn kernel_polys<F: Field>(f: &F, e: &Curve<F::E>, ell: u64, rng: &mut Rng) -> Vec<Poly<F>> {
    let psi = poly::monic(f, &division_poly(f, e, ell as usize));
    let facs = poly::factor_squarefree(f, &psi, rng);
    let d = ((ell - 1) / 2) as usize;
    let mut out = vec![];
    let degs: Vec<usize> = facs.iter().map(|p| p.len() - 1).collect();
    // enumerate subsets with total degree d
    fn rec(i: usize, left: usize, degs: &[usize], cur: &mut Vec<usize>, res: &mut Vec<Vec<usize>>) {
        if left == 0 {
            res.push(cur.clone());
            return;
        }
        if i == degs.len() {
            return;
        }
        if degs[i] <= left {
            cur.push(i);
            rec(i + 1, left - degs[i], degs, cur, res);
            cur.pop();
        }
        rec(i + 1, left, degs, cur, res);
    }
    let mut subsets = vec![];
    rec(0, d, &degs, &mut vec![], &mut subsets);
    for s in subsets {
        let mut h = vec![f.one()];
        for &i in &s {
            h = poly::mul(f, &h, &facs[i]);
        }
        let iso = kohel(f, e, &h, ell);
        if check_x_identity(f, &iso, rng, 4) {
            out.push(h);
        }
    }
    out
}
