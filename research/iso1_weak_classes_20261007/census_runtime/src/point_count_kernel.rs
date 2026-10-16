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

// ── Certifying counter (added 2026-10-16; not part of the frozen kernel) ──
//
// `#E(F_Q)` lies in the Hasse interval `[Q+1-2p³, Q+1+2p³]` (here `Q = p⁶`,
// so `√Q = p³` exactly) and is divisible by 4 because all three 2-torsion
// points are rational; the same holds for the quadratic twist, whose order
// is `2Q+2-N`. The counter keeps the full set of candidates `N` compatible
// with every point tested, on the curve and on its twist, and returns `N`
// only when exactly one candidate survives. A returned order is therefore
// proved, not sampled: the true order always survives every filter.

/// Outcome of [`certified_order`]: the unique surviving order and the
/// number of curve and twist points used to isolate it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Certified {
    pub order: u128,
    pub curve_points: u32,
    pub twist_points: u32,
}

/// Every `k` in `[klo, khi]` with `k·r = O`, or `None` when there are more
/// than `cap` of them (an uninformative point). Scans the whole range.
fn all_annihilators(ec: &EllE, r: &PtE6, klo: u128, khi: u128, cap: usize) -> Option<Vec<u128>> {
    if r.inf {
        return None;
    }
    let width = khi - klo;
    let m = isqrt_u128(width) + 1;
    let mut table: HashMap<PtE6, u128> = HashMap::with_capacity(m as usize);
    let mut jp = *r;
    for j in 1..m {
        if jp.inf {
            // ord(r) = j < m: list its multiples directly.
            let first = klo.div_ceil(j) * j;
            let count = if first > khi { 0 } else { (khi - first) / j + 1 };
            if count as usize > cap {
                return None;
            }
            return Some((0..count).map(|i| first + i * j).collect());
        }
        table.entry(jp).or_insert(j);
        jp = ec.add(&jp, r);
    }
    // ord(r) >= m, so each giant block of m consecutive k holds at most one
    // annihilator, found either as T = O or as T = -j·r.
    let giant = ec.mul_u128(r, m);
    let mut t = ec.mul_u128(r, klo);
    let mut out = Vec::new();
    let mut base = klo;
    while base <= khi {
        let hit = if t.inf {
            Some(base)
        } else {
            table.get(&ec.neg(&t)).map(|&j| base + j)
        };
        if let Some(k) = hit {
            if k <= khi {
                out.push(k);
                if out.len() > cap {
                    return None;
                }
            }
        }
        t = ec.add(&t, &giant);
        base += m;
    }
    Some(out)
}

/// The order of `y² = x³ + a₂x² + a₄x` over `F_{p⁶}`, certified unique among
/// Hasse-interval multiples of 4 compatible with the tested points; `None`
/// if `max_points` points on each side leave more than one candidate.
pub fn certified_order(
    f: &Fq3,
    a2: E6,
    a4: E6,
    nonsquare: &E6,
    rng: &mut StdRng,
    max_points: u32,
) -> Option<Certified> {
    let p = f.f.p as u128;
    let q3 = p.pow(6);
    let s = p.pow(3);
    let (lo, hi) = (q3 + 1 - 2 * s, q3 + 1 + 2 * s);
    let (klo, khi) = (lo.div_ceil(4), hi / 4);
    let ec = EllE::from_a2_a4(f, a2, a4);
    let d = nonsquare;
    let tw = EllE::from_a2_a4(f, f.mul(d, &a2), f.mul(&f.sq(d), &a4));
    let mut cand: Option<Vec<u128>> = None;
    let (mut cp, mut tp) = (0_u32, 0_u32);
    let done = |c: &Option<Vec<u128>>| matches!(c, Some(v) if v.len() == 1);
    while cp < max_points || tp < max_points {
        // Alternate curve and twist points; a curve point is tried first.
        let use_twist = (tp < cp || cp >= max_points) && tp < max_points;
        let (curve, n_of) = if use_twist { (&tw, true) } else { (&ec, false) };
        let pt = curve.random_point(rng);
        if use_twist {
            tp += 1;
        } else {
            cp += 1;
        }
        // Orders of the curve: N; of the twist: 2Q+2-N. Both are 0 mod 4.
        let to_n = |m: u128| if n_of { 2 * q3 + 2 - m } else { m };
        cand = match cand {
            Some(v) => Some(
                v.into_iter()
                    .filter(|&n| curve.mul_u128(&pt, to_n(n)).inf)
                    .collect(),
            ),
            None => {
                let r = curve.add(&curve.add(&pt, &pt), &curve.add(&pt, &pt));
                all_annihilators(curve, &r, klo, khi, 4096)
                    .map(|ks| ks.into_iter().map(|k| to_n(4 * k)).collect())
            }
        };
        if let Some(v) = &cand {
            assert!(!v.is_empty(), "true order eliminated: arithmetic error");
        }
        if done(&cand) {
            let order = cand.unwrap()[0];
            return Some(Certified {
                order,
                curve_points: cp,
                twist_points: tp,
            });
        }
    }
    None
}
