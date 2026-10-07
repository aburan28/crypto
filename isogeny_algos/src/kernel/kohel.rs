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
        f.sub(
            f.mul(e1, f.mul(e1, e1)),
            f.mul(f.from_u64(3), f.mul(e1, e2)),
        ),
        f.mul(f.from_u64(3), e3),
    );
    (p1, p2, p3)
}

pub fn codomain_from_sums<F: Field>(
    f: &F,
    e: &Curve<F::E>,
    n: usize,
    p: (F::E, F::E, F::E),
) -> Curve<F::E> {
    let nn = f.from_u64(n as u64);
    let t = f.add(
        f.mul(f.from_u64(6), p.1),
        f.mul(f.from_u64(2), f.mul(nn, e.a)),
    );
    let w = f.add(
        f.add(
            f.mul(f.from_u64(10), p.2),
            f.mul(f.from_u64(6), f.mul(e.a, p.0)),
        ),
        f.mul(f.from_u64(4), f.mul(nn, e.b)),
    );
    Curve::new(
        f.sub(e.a, f.mul(f.from_u64(5), t)),
        f.sub(e.b, f.mul(f.from_u64(7), w)),
    )
}

/// Kohel's isogeny for the kernel polynomial `h` = prod over S of (x - x_Q), where S holds one
/// representative per +-pair of non-zero kernel points (2-torsion points once). Handles odd and
/// even degree: the 2-torsion part h2 = gcd(h, x^3+Ax+B) is split off and enters with simple
/// poles, the rest hr with double poles; deg = 1 + deg h2 + 2 deg hr.
pub fn kohel<F: Field>(f: &F, e: &Curve<F::E>, h: &Poly<F>, ell: u64) -> RatIsogeny<F> {
    let c = |v: u64| f.from_u64(v);
    let fx: Poly<F> = vec![e.b, e.a, f.zero(), f.one()];
    let h2 = poly::gcd(f, h, &fx);
    let hr = poly::divrem(f, h, &h2).0;
    let (k, r) = (h2.len() - 1, hr.len() - 1);
    debug_assert_eq!(
        ell,
        1 + k as u64 + 2 * r as u64,
        "kernel polynomial degree does not match the isogeny degree"
    );
    let (p2s, pr) = (power_sums(f, &h2), power_sums(f, &hr));
    let (kk, rr) = (c(k as u64), c(r as u64));
    // t = sum_{G2}(3x^2+A) + sum_R(6x^2+2A);  w = sum_{G2}(3x^3+Ax) + sum_R(10x^3+6Ax+4B)
    let t = f.add(
        f.add(f.mul(c(3), p2s.1), f.mul(kk, e.a)),
        f.add(f.mul(c(6), pr.1), f.mul(f.mul(c(2), rr), e.a)),
    );
    let w = f.add(
        f.add(
            f.add(f.mul(c(3), p2s.2), f.mul(e.a, p2s.0)),
            f.mul(c(10), pr.2),
        ),
        f.add(f.mul(c(6), f.mul(e.a, pr.0)), f.mul(c(4), f.mul(rr, e.b))),
    );
    let cod = Curve::new(f.sub(e.a, f.mul(c(5), t)), f.sub(e.b, f.mul(c(7), w)));
    // f(x) = x + f_R + f_2 with
    //   f_R = -(6x^2+2A) s1_R + 4F s2_R + 2(r x - p1_R),  s1_R = hr'/hr, s2_R = (hr'^2 - hr hr'')/hr^2
    //   f_2 = F'(x) h2'/h2 - 3(k x + p1_2)
    // N = f * hr^2 h2.
    let hr1 = poly::derivative(f, &hr);
    let hr2 = poly::derivative(f, &hr1);
    let h2d = poly::derivative(f, &h2);
    let hr_sq = poly::mul(f, &hr, &hr);
    let den = poly::mul(f, &hr_sq, &h2);
    let q = poly::sub(f, &poly::mul(f, &hr1, &hr1), &poly::mul(f, &hr, &hr2));
    let g1 = vec![f.mul(c(2), e.a), f.zero(), c(6)]; // 6x^2 + 2A
    let fp1 = vec![e.a, f.zero(), c(3)]; // F'(x) = 3x^2 + A
    let lin_r = vec![f.neg(f.mul(c(2), pr.0)), f.mul(c(2), rr)]; // 2(r x - p1_R)
    let lin_2 = vec![f.neg(f.mul(c(3), p2s.0)), f.neg(f.mul(c(3), kk))]; // -3(k x + p1_2)
    let xx = vec![f.zero(), f.one()];
    let mut num = poly::mul(f, &xx, &den);
    num = poly::sub(
        f,
        &num,
        &poly::mul(f, &g1, &poly::mul(f, &poly::mul(f, &hr, &hr1), &h2)),
    );
    num = poly::add(
        f,
        &num,
        &poly::scale(f, &poly::mul(f, &poly::mul(f, &fx, &q), &h2), c(4)),
    );
    num = poly::add(f, &num, &poly::mul(f, &lin_r, &den));
    num = poly::add(f, &num, &poly::mul(f, &fp1, &poly::mul(f, &hr_sq, &h2d)));
    num = poly::add(f, &num, &poly::mul(f, &lin_2, &den));
    RatIsogeny {
        dom: *e,
        cod,
        deg: ell,
        ker: h.clone(),
        num,
        den,
    }
}
