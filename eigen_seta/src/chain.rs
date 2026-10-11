//! Decapsulation chain with a fixed public schedule, and the key-generation
//! chain.  See RESEARCH note §4: the eigen-structure forces
//!   graph, prodker, glue, rich x (e-6), split, prodker, graph
//! and key generation selects keys for which this schedule holds, so the
//! sequence of operations never depends on the secret key or the ciphertext.

use crate::ct::{Choice, ConditionallySelectable};
use crate::field::Field;
use crate::kani::{div_dbl_iter, div_dbl_iter_batch, Div, Glue, Richelot, Split};
use crate::strategy::optimal;
use crate::weier::{Isog2, WCurve, WIso, WPoint};

pub type Pair<F> = (WPoint<F>, WPoint<F>);

pub fn graph_identity<F: Field>(e: &WCurve<F>, p: &Pair<F>) -> Pair<F> {
    (e.add(&p.0, &p.1), e.sub(&p.0, &p.1))
}

/// Product-kernel step: each factor gets one 2-torsion point of the kernel,
/// selected in constant time (exactly one of the two generators is nonzero
/// on each factor for honest inputs).
pub fn prodker_step<F: Field>(
    e1: &WCurve<F>,
    e2: &WCurve<F>,
    k1: &Pair<F>,
    k2: &Pair<F>,
    h: usize,
) -> (Isog2<F>, Isog2<F>) {
    let a1 = e1.dbl_iter(&k1.0, h);
    let a2 = e1.dbl_iter(&k2.0, h);
    let b1 = e2.dbl_iter(&k1.1, h);
    let b2 = e2.dbl_iter(&k2.1, h);
    let kap1 = WPoint::conditional_select(&a1, &a2, a1.inf);
    let kap2 = WPoint::conditional_select(&b1, &b2, b1.inf);
    (Isog2::new(e1, &kap1.x), Isog2::new(e2, &kap2.x))
}

pub fn map_pair<F: Field>(f: &(Isog2<F>, Isog2<F>), p: &Pair<F>) -> Pair<F> {
    (f.0.eval(&p.0), f.1.eval(&p.1))
}

/// x-coordinates (x(P), x(Q), x(P+Q)) of a 2-torsion basis, the third from
/// the root sum of x^3 + a2 x^2 + ...
pub fn two_torsion_xs<F: Field>(e: &WCurve<F>, p: &WPoint<F>, q: &WPoint<F>) -> [F; 3] {
    [p.x, q.x, -e.a2 - p.x - q.x]
}

/// Jacobian phase: `r` Richelot steps along the kernel <D1, D2> (order >= 2^r),
/// traversed with an optimal strategy; the generators and `extra` divisors
/// are pushed through every step.  Returns the final h.
pub fn jacobian_phase<F: Field>(
    h0: [F; 7],
    gens: &mut [Div<F>; 2],
    order_log: usize,
    r: usize,
    extra: &mut [Div<F>],
) -> [F; 7] {
    let mut h = h0;
    if r == 0 {
        return h;
    }
    let strategy = optimal(r, 0.6, 1.0);
    // leaf generator of order 2^r
    let g0v = div_dbl_iter_batch(&h, &gens[..], order_log - r);
    let mut stack: Vec<([Div<F>; 2], usize)> = vec![([g0v[0], g0v[1]], 0)];
    let mut idx = 0;
    for row in 0..r {
        loop {
            let (q, k) = *stack.last().unwrap();
            if k >= r - 1 - row {
                break;
            }
            let s = strategy[idx];
            idx += 1;
            let m = div_dbl_iter_batch(&h, &q, s);
            stack.push(([m[0], m[1]], k + s));
        }
        let (leaf, _) = stack.pop().unwrap();
        let phi = Richelot::new(&h, &leaf[0], &leaf[1]);
        // push everything through in one batch
        let mut all: Vec<Div<F>> = Vec::with_capacity(2 * stack.len() + 2 + extra.len());
        for e in stack.iter() {
            all.push(e.0[0]);
            all.push(e.0[1]);
        }
        all.push(gens[0]);
        all.push(gens[1]);
        all.extend_from_slice(extra);
        let mapped = phi.map_batch(&all);
        let mut it = mapped.into_iter();
        for e in stack.iter_mut() {
            e.0 = [it.next().unwrap(), it.next().unwrap()];
        }
        gens[0] = it.next().unwrap();
        gens[1] = it.next().unwrap();
        for d in extra.iter_mut() {
            *d = it.next().unwrap();
        }
        h = phi.hnew;
    }
    h
}

/// The decapsulation (2^e, 2^e)-chain on E x E with kernel <K1, K2> and the
/// fixed schedule.  Returns the codomain curve (both factors are the same
/// curve after the final graph step) and the images of `payload`.
pub fn decaps_chain<F: Field>(ec: &WCurve<F>, k1: &Pair<F>, k2: &Pair<F>, e: usize, payload: &[Pair<F>]) -> (WCurve<F>, Vec<Pair<F>>) {
    assert!(e >= 8);
    // step 1: graph of the identity  (x, y) -> (x + y, x - y)
    let mut k1 = graph_identity(ec, k1);
    let mut k2 = graph_identity(ec, k2);
    let mut pay: Vec<Pair<F>> = payload.iter().map(|p| graph_identity(ec, p)).collect();
    // step 2: product kernel (kernel now of order 2^(e-1))
    let f = prodker_step(ec, ec, &k1, &k2, e - 2);
    k1 = map_pair(&f, &k1);
    k2 = map_pair(&f, &k2);
    pay = pay.iter().map(|p| map_pair(&f, p)).collect();
    let (c1, c2) = (f.0.codomain, f.1.codomain);
    // step 3: gluing (kernel of order 2^(e-2))
    let a1 = c1.dbl_iter(&k1.0, e - 3);
    let a2 = c1.dbl_iter(&k2.0, e - 3);
    let b1 = c2.dbl_iter(&k1.1, e - 3);
    let b2 = c2.dbl_iter(&k2.1, e - 3);
    let g = Glue::new(two_torsion_xs(&c1, &a1, &a2), two_torsion_xs(&c2, &b1, &b2));
    let mut gens = [g.map(&k1.0, &k1.1), g.map(&k2.0, &k2.1)];
    let mut extra: Vec<Div<F>> = pay.iter().map(|p| g.map(&p.0, &p.1)).collect();
    // steps 4 .. e-3: Richelot (kernel of order 2^(e-3) after gluing)
    let r = e - 6;
    let h = jacobian_phase(g.h, &mut gens, e - 3, r, &mut extra);
    // step e-2: splitting (kernel of order 2^3 now; its 2-torsion is 4*gens)
    let s = Split::new(&h, &div_dbl_iter(&h, &gens[0], 2), &div_dbl_iter(&h, &gens[1], 2));
    let mut k1 = s.map(&gens[0]);
    let mut k2 = s.map(&gens[1]);
    let mut pay: Vec<Pair<F>> = extra.iter().map(|d| s.map(d)).collect();
    let (e3, e4) = (s.e1, s.e2);
    // step e-1: product kernel (kernel of order 4)
    let f = prodker_step(&e3, &e4, &k1, &k2, 1);
    k1 = map_pair(&f, &k1);
    k2 = map_pair(&f, &k2);
    let _ = k2;
    pay = pay.iter().map(|p| map_pair(&f, p)).collect();
    let (e5, e6) = (f.0.codomain, f.1.codomain);
    // step e: graph of an isomorphism e6 -> e5 on the 2-torsion
    let tau = WIso::find(&e6, &e5);
    let _ = k1;
    let out: Vec<Pair<F>> = pay
        .iter()
        .map(|p| {
            let t = tau.eval(&p.1);
            (e5.add(&p.0, &t), e5.sub(&p.0, &t))
        })
        .collect();
    (e5, out)
}

/// Key-generation chain (variable time is acceptable here): schedule
/// glue, rich x (e-2), split, starting from E x E with kernel <K1, K2> of
/// order 2^e.  Returns the two codomain curves and the payload images.
pub fn keygen_chain<F: Field>(e0: &WCurve<F>, k1: &Pair<F>, k2: &Pair<F>, e: usize, payload: &[Pair<F>]) -> (WCurve<F>, WCurve<F>, Vec<Pair<F>>) {
    let a1 = e0.dbl_iter(&k1.0, e - 1);
    let a2 = e0.dbl_iter(&k2.0, e - 1);
    let b1 = e0.dbl_iter(&k1.1, e - 1);
    let b2 = e0.dbl_iter(&k2.1, e - 1);
    let g = Glue::new(two_torsion_xs(e0, &a1, &a2), two_torsion_xs(e0, &b1, &b2));
    let mut gens = [g.map(&k1.0, &k1.1), g.map(&k2.0, &k2.1)];
    let mut extra: Vec<Div<F>> = payload.iter().map(|p| g.map(&p.0, &p.1)).collect();
    let h = jacobian_phase(g.h, &mut gens, e - 1, e - 2, &mut extra);
    let s = Split::new(&h, &gens[0], &gens[1]);
    let pay: Vec<Pair<F>> = extra.iter().map(|d| s.map(d)).collect();
    (s.e1, s.e2, pay)
}

/// Does the first (2,2)-step of a keygen kernel glue (no product/graph
/// structure)?  Variable time; used to reject keygen candidates.
pub fn first_step_glues<F: Field>(e0: &WCurve<F>, k1: &Pair<F>, k2: &Pair<F>, e: usize) -> bool {
    let a1 = e0.dbl_iter(&k1.0, e - 1);
    let a2 = e0.dbl_iter(&k2.0, e - 1);
    let b1 = e0.dbl_iter(&k1.1, e - 1);
    let b2 = e0.dbl_iter(&k2.1, e - 1);
    let fin = |p: &WPoint<F>| !bool::from(p.inf);
    // graph of the identity on E[2] (alpha = 1 on E[2]) gives a product step
    fin(&a1) && fin(&a2) && fin(&b1) && fin(&b2) && !(bool::from(a1.x.ct_eq(&b1.x)) && bool::from(a2.x.ct_eq(&b2.x)))
}

pub(crate) fn _c(_: Choice) {}
