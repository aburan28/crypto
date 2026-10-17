//! Rational 3-isogenies of y^2 = x^3 + A x + B over GF(p), p > 3, and the
//! 3-volcano position of a curve. psi_3 = 3x^4 + 6A x^2 + 12B x - A^2; each
//! rational root x0 is a kernel {O, P, -P}. Velu (one representative P):
//! v = 6 x0^2 + 2A, w = 4 (x0^3 + A x0 + B) + x0 v, A' = A - 5v, B' = B - 7w.
//! Vertices are j-invariants (j != 0, 1728 is checked): with the trace
//! fixed, j determines the F_p-isomorphism class.

use crate::fp::Fp;
use std::collections::{HashMap, VecDeque};

type P = Vec<u64>;

fn trim(mut a: P) -> P {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
    a
}
fn deg(a: &P) -> isize {
    if a.len() == 1 && a[0] == 0 {
        -1
    } else {
        a.len() as isize - 1
    }
}
fn sub(k: &Fp, a: &P, b: &P) -> P {
    let n = a.len().max(b.len());
    let mut r = vec![0u64; n];
    for i in 0..n {
        let x = *a.get(i).unwrap_or(&0);
        let y = *b.get(i).unwrap_or(&0);
        r[i] = k.sub(x, y);
    }
    trim(r)
}
fn mul(k: &Fp, a: &P, b: &P) -> P {
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = k.add(r[i + j], k.mul(x, y));
        }
    }
    trim(r)
}
fn rem(k: &Fp, a: &P, m: &P) -> P {
    let dm = deg(m) as usize;
    let inv = k.inv(m[dm]);
    let mut r = a.clone();
    while deg(&r) >= dm as isize {
        let dr = deg(&r) as usize;
        let c = k.mul(r[dr], inv);
        let s = dr - dm;
        for (i, &mi) in m.iter().enumerate() {
            r[s + i] = k.sub(r[s + i], k.mul(c, mi));
        }
        r = trim(r);
    }
    r
}
fn monic(k: &Fp, a: P) -> P {
    let d = deg(&a);
    if d < 0 {
        return a;
    }
    let inv = k.inv(a[d as usize]);
    a.iter().map(|&c| k.mul(c, inv)).collect()
}
fn gcd(k: &Fp, a: &P, b: &P) -> P {
    let (mut a, mut b) = (a.clone(), b.clone());
    while deg(&b) >= 0 {
        let r = rem(k, &a, &b);
        a = b;
        b = r;
    }
    monic(k, a)
}
fn powmod(k: &Fp, base: &P, mut e: u128, m: &P) -> P {
    let mut r: P = vec![1];
    let mut b = rem(k, base, m);
    while e > 0 {
        if e & 1 == 1 {
            r = rem(k, &mul(k, &r, &b), m);
        }
        b = rem(k, &mul(k, &b, &b), m);
        e >>= 1;
    }
    r
}
fn div_exact(k: &Fp, a: &P, b: &P) -> P {
    let db = deg(b) as usize;
    let inv = k.inv(b[db]);
    let mut r = a.clone();
    let mut q = vec![0u64; deg(a) as usize - db + 1];
    while deg(&r) >= db as isize {
        let dr = deg(&r) as usize;
        let c = k.mul(r[dr], inv);
        q[dr - db] = c;
        for (i, &bi) in b.iter().enumerate() {
            r[dr - db + i] = k.sub(r[dr - db + i], k.mul(c, bi));
        }
        r = trim(r);
    }
    assert!(deg(&r) < 0, "inexact division");
    trim(q)
}

/// Distinct roots in F_p of f (Cantor-Zassenhaus).
pub fn roots(k: &Fp, f: &P) -> Vec<u64> {
    let f = monic(k, trim(f.clone()));
    let xp = powmod(k, &vec![0, 1], k.p as u128, &f);
    let g = gcd(k, &f, &sub(k, &xp, &vec![0, 1]));
    let mut out = Vec::new();
    split(k, &g, &mut out, 1);
    out.sort_unstable();
    out
}

fn split(k: &Fp, g: &P, out: &mut Vec<u64>, mut delta: u64) {
    let d = deg(g);
    if d <= 0 {
        return;
    }
    if d == 1 {
        out.push(k.sub(0, g[0]));
        return;
    }
    loop {
        delta += 1;
        let h = powmod(k, &vec![delta % k.p, 1], (k.p as u128 - 1) / 2, g);
        let h = gcd(k, g, &sub(k, &h, &vec![1]));
        let dh = deg(&h);
        if dh > 0 && dh < d {
            let q = div_exact(k, g, &h);
            split(k, &h, out, delta * 7919);
            split(k, &q, out, delta * 104729);
            return;
        }
    }
}

pub fn j_invariant(k: &Fp, a: u64, b: u64) -> u64 {
    let a3 = k.mul(4, k.mul(a, k.mul(a, a)));
    let den = k.add(a3, k.mul(27, k.mul(b, b)));
    k.mul(k.mul(1728, a3), k.inv(den))
}

/// The rational 3-isogenous curves (A', B').
pub fn neighbours(k: &Fp, a: u64, b: u64) -> Vec<(u64, u64)> {
    let psi3: P = vec![k.sub(0, k.mul(a, a)), k.mul(12, b), k.mul(6, a), 0, 3];
    roots(k, &psi3)
        .into_iter()
        .map(|x0| {
            let v = k.add(k.mul(6, k.mul(x0, x0)), k.mul(2, a));
            let fx = k.add(k.mul(k.add(k.mul(x0, x0), a), x0), b);
            let w = k.add(k.mul(4, fx), k.mul(x0, v));
            (k.sub(a, k.mul(5, v)), k.sub(b, k.mul(7, w)))
        })
        .collect()
}

#[derive(Debug)]
pub struct Volcano {
    pub height: u32,
    pub root_neighbours: usize,
    pub distance_to_floor: Option<u32>,
    pub level_from_crater: Option<u32>,
    pub vertices_explored: usize,
}

/// BFS over the 3-isogeny graph to the nearest floor vertex (one rational
/// 3-isogeny); the distance is the root's height above the floor.
pub fn volcano3(k: &Fp, a: u64, b: u64, height: u32) -> Volcano {
    let j0 = j_invariant(k, a, b);
    assert!(j0 != 0 && j0 != 1728 % k.p);
    let n0 = neighbours(k, a, b);
    let mut seen: HashMap<u64, u32> = HashMap::new();
    seen.insert(j0, 0);
    let mut q = VecDeque::new();
    q.push_back((j0, n0.clone()));
    let mut found = None;
    while let Some((jv, nv)) = q.pop_front() {
        let dv = seen[&jv];
        if height > 0 && nv.len() == 1 {
            found = Some(dv);
            break;
        }
        if height == 0 || dv >= height {
            continue;
        }
        for (a2, b2) in nv {
            let jw = j_invariant(k, a2, b2);
            assert!(jw != 0 && jw != 1728 % k.p);
            if !seen.contains_key(&jw) {
                seen.insert(jw, dv + 1);
                q.push_back((jw, neighbours(k, a2, b2)));
            }
        }
    }
    if height == 0 {
        found = Some(0);
    }
    Volcano {
        height,
        root_neighbours: n0.len(),
        distance_to_floor: found,
        level_from_crater: found.map(|d| height - d),
        vertices_explored: seen.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::PCurve;

    fn naive(k: &Fp, a: u64, b: u64) -> u64 {
        let mut n = 1;
        for x in 0..k.p {
            let r = k.add(k.mul(k.add(k.mul(x, x), a), x), b);
            n += if r == 0 { 1 } else if k.is_square(r) { 2 } else { 0 };
        }
        n
    }

    #[test]
    fn neighbours_have_equal_order() {
        let k = Fp { p: 10_007 };
        let mut tested = 0;
        for b in 1..60u64 {
            let a = k.p - 3;
            let _ = PCurve { k: k.clone(), a, b };
            let n = naive(&k, a, b);
            for (a2, b2) in neighbours(&k, a, b) {
                assert_eq!(naive(&k, a2, b2), n);
                tested += 1;
            }
        }
        assert!(tested > 10);
    }
}
