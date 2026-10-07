//! Kohel (thesis, 1996) formulas: codomain and rational x-map from the kernel
//! polynomial h(x) = prod (x - x_Q) alone; no kernel points needed.
use crate::curve::*;
use crate::field::Field;
use crate::poly::{self, Poly};

/// Power sums p1,p2,p3 of the roots from the top coefficients of monic h of degree n.
pub fn power_sums<F: Field>(f: &F, h: &Poly<F>) -> (F::E, F::E, F::E) {
    let n = h.len() - 1;
    let c = |k: usize| if k <= n { h[n - k] } else { f.zero() }; // coefficient of x^{n-k}
    // h = x^n - e1 x^{n-1} + e2 x^{n-2} - e3 x^{n-3}
    let e1 = f.neg(c(1));
    let e2 = c(2);
    let e3 = f.neg(c(3));
    let p1 = e1;
    let p2 = f.sub(f.mul(e1, e1), f.mul(f.from_u64(2), e2));
    let p3 = f.add(
        f.sub(f.mul(e1, f.mul(e1, e1)), f.mul(f.from_u64(3), f.mul(e1, e2))),
        f.mul(f.from_u64(3), e3),
    );
    (p1, p2, p3)
}

pub fn codomain_from_sums<F: Field>(f: &F, e: &Curve<F::E>, n: usize, p: (F::E, F::E, F::E)) -> Curve<F::E> {
    let nn = f.from_u64(n as u64);
    let t = f.add(f.mul(f.from_u64(6), p.1), f.mul(f.from_u64(2), f.mul(nn, e.a)));
    let w = f.add(
        f.add(f.mul(f.from_u64(10), p.2), f.mul(f.from_u64(6), f.mul(e.a, p.0))),
        f.mul(f.from_u64(4), f.mul(nn, e.b)),
    );
    Curve::new(f.sub(e.a, f.mul(f.from_u64(5), t)), f.sub(e.b, f.mul(f.from_u64(7), w)))
}

/// Kohel's isogeny for kernel polynomial `h` (monic, degree (l-1)/2).
pub fn kohel<F: Field>(f: &F, e: &Curve<F::E>, h: &Poly<F>, ell: u64) -> RatIsogeny<F> {
    let n = h.len() - 1;
    let ps = power_sums(f, h);
    let cod = codomain_from_sums(f, e, n, ps);
    // f(x) = x - (6x^2+2A) s1 + 4F(x) s2 + 2(n x - p1), s1 = h'/h, s2 = (h'^2 - h h'')/h^2
    // => N = [x + 2(n x - p1)] h^2 - (6x^2+2A) h h' + 4F (h'^2 - h h'')
    let c = |v: u64| f.from_u64(v);
    let h1 = poly::derivative(f, h);
    let h2 = poly::derivative(f, &h1);
    let hh = poly::mul(f, h, h);
    let hh1 = poly::mul(f, h, &h1);
    let q = poly::sub(f, &poly::mul(f, &h1, &h1), &poly::mul(f, h, &h2));
    let nn = c(n as u64);
    let lin = vec![f.neg(f.mul(c(2), ps.0)), f.add(c(1), f.mul(c(2), nn))];
    let quad = vec![f.mul(c(2), e.a), f.zero(), c(6)];
    let fx = vec![e.b, e.a, f.zero(), f.one()];
    let mut num = poly::mul(f, &lin, &hh);
    num = poly::sub(f, &num, &poly::mul(f, &quad, &hh1));
    num = poly::add(f, &num, &poly::scale(f, &poly::mul(f, &fx, &q), c(4)));
    RatIsogeny { dom: *e, cod, deg: ell, ker: h.clone(), num, den: hh }
}
