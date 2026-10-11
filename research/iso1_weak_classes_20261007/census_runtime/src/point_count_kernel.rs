//! Randomized two-point Hasse-interval counter frozen from 62becf9572fe74cbe8b3d8cebee3bf8a240708a1.
use crate::field_kernel::{EllE, Fld, Fq3, PtE6, E6};
use rand::rngs::StdRng;
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve2 {
    pub e: [E6; 3],
}
pub fn norm_q3(f: &Fq3, z: &E6) -> E6 {
    f.mul(z, &f.mul(&f.sigma(z), &f.sigma(&f.sigma(z))))
}
impl Curve2 {
    pub fn weak_by_norms(&self, f: &Fq3) -> bool {
        let n01 = norm_q3(f, &f.sub(&self.e[1], &self.e[0]));
        let n02 = norm_q3(f, &f.sub(&self.e[2], &self.e[0]));
        let n12 = norm_q3(f, &f.sub(&self.e[2], &self.e[1]));
        // (e₀; e₁, e₂): N(e₂−e₀) = N(e₁−e₀); (e₁; e₂, e₀): N(e₀−e₁) = N(e₂−e₁);
        // (e₂; e₀, e₁): N(e₁−e₂) = N(e₀−e₂)
        n02 == n01 || n12 == f.neg(&n01) || n12 == n02
    }
}
pub fn random_curve(f: &Fq3, rng: &mut StdRng) -> Curve2 {
    loop {
        let e = [f.random(rng), f.random(rng), f.random(rng)];
        if e[0] != e[1] && e[1] != e[2] && e[0] != e[2] {
            return Curve2 { e };
        }
    }
}
fn isqrt_u128(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u128;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}
pub fn curve_order(f: &Fq3, c: &Curve2, rng: &mut StdRng) -> u128 {
    let u = f.sub(&c.e[1], &c.e[0]);
    let v = f.sub(&c.e[2], &c.e[0]);
    let ec = EllE::from_a2_a4(f, f.neg(&f.add(&u, &v)), f.mul(&u, &v));
    ell_order(f, &ec, rng)
}
pub fn ell_order(f: &Fq3, ec: &EllE, rng: &mut StdRng) -> u128 {
    let p = f.f.p as u128;
    let q3 = p.pow(6);
    let two_sqrt = 2 * isqrt_u128(q3) + 2;
    let lo = q3 + 1 - two_sqrt;
    let width = 2 * two_sqrt;
    let steps = isqrt_u128(width) + 1;
    let order_of = |pt: &PtE6| -> Option<u128> {
        let mut table: HashMap<PtE6, u128> = HashMap::new();
        let mut jp = PtE6::INF;
        for j in 0..steps {
            if !jp.inf {
                table.entry(jp).or_insert(j);
            }
            jp = ec.add(&jp, pt);
        }
        let giant = ec.mul_u128(pt, steps);
        let mut t = ec.mul_u128(pt, lo);
        let mut i = 0u128;
        while i * steps <= width + steps {
            if t.inf {
                return Some(lo + i * steps);
            }
            if let Some(&j) = table.get(&ec.neg(&t)) {
                return Some(lo + i * steps + j);
            }
            t = ec.add(&t, &giant);
            i += 1;
        }
        None
    };
    loop {
        let p1 = ec.random_point(rng);
        let Some(m) = order_of(&p1) else { continue };
        let p2 = ec.random_point(rng);
        if ec.mul_u128(&p2, m).inf {
            return m;
        }
    }
}
