//! Rational 3-isogenies of ordinary binary curves E_b: y^2 + xy = x^3 + x^2 + b
//! over GF(2^m), m odd, and their 3-volcano position.
//!
//! The 3-division polynomial of y^2 + xy = x^3 + a2 x^2 + a6 is
//! psi_3 = x^4 + x^3 + a6. Each rational root x0 is the x-coordinate of a
//! kernel {O, P, -P}; Velu's formulas in characteristic 2 (a1 = 1, a3 = a4 = 0)
//! give A4 = x0, A6 = a6 + x0, and the substitution y -> y + A4 returns the
//! model y^2 + xy = x^3 + a2 x^2 + (a6 + x0 + x0^2). The isogenous curve has
//! the same order, so it keeps a2 = 1 (Tr(a2) = 1 for odd m).

use crate::f2m::F2m;
use std::collections::{HashMap, VecDeque};

pub type Poly = Vec<u64>; // coefficients, constant term first

pub fn trim(mut a: Poly) -> Poly {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
    if a.is_empty() {
        a.push(0);
    }
    a
}

pub fn deg(a: &Poly) -> isize {
    if a.len() == 1 && a[0] == 0 {
        -1
    } else {
        a.len() as isize - 1
    }
}

pub fn add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0u64; n];
    for (i, &x) in a.iter().enumerate() {
        r[i] ^= x;
    }
    for (i, &x) in b.iter().enumerate() {
        r[i] ^= x;
    }
    trim(r)
}

pub fn mul(k: &F2m, a: &Poly, b: &Poly) -> Poly {
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] ^= k.mul(x, y);
        }
    }
    trim(r)
}

pub fn divrem(k: &F2m, a: &Poly, b: &Poly) -> (Poly, Poly) {
    let db = deg(b);
    assert!(db >= 0);
    let inv_lead = k.inv(b[db as usize]);
    let mut r = a.clone();
    let mut q = vec![0u64; (deg(a) - db + 1).max(1) as usize];
    while deg(&r) >= db {
        let dr = deg(&r) as usize;
        let c = k.mul(r[dr], inv_lead);
        let s = dr - db as usize;
        q[s] ^= c;
        for (i, &bi) in b.iter().enumerate() {
            r[s + i] ^= k.mul(c, bi);
        }
        r = trim(r);
    }
    (trim(q), r)
}

pub fn gcd(k: &F2m, a: &Poly, b: &Poly) -> Poly {
    let (mut a, mut b) = (a.clone(), b.clone());
    while deg(&b) >= 0 {
        let (_, r) = divrem(k, &a, &b);
        a = b;
        b = r;
    }
    // monic
    let d = deg(&a);
    if d >= 0 {
        let inv = k.inv(a[d as usize]);
        a = a.iter().map(|&c| k.mul(c, inv)).collect();
    }
    a
}

pub fn mulmod(k: &F2m, a: &Poly, b: &Poly, m: &Poly) -> Poly {
    divrem(k, &mul(k, a, b), m).1
}

/// Distinct roots in GF(2^m) of a squarefree polynomial.
pub fn roots(k: &F2m, f: &Poly) -> Vec<u64> {
    let f = trim(f.clone());
    // g = gcd(f, x^(2^m) - x)
    let x: Poly = vec![0, 1];
    let mut xp = x.clone();
    for _ in 0..k.m {
        xp = mulmod(k, &xp, &xp, &f);
    }
    let g = gcd(k, &f, &add(&xp, &x));
    let mut out = Vec::new();
    split(k, &g, &mut out, 1);
    out.sort_unstable();
    out
}

fn split(k: &F2m, g: &Poly, out: &mut Vec<u64>, mut seed: u64) {
    let d = deg(g);
    if d <= 0 {
        return;
    }
    if d == 1 {
        // monic: x + c
        let inv = k.inv(g[1]);
        out.push(k.mul(g[0], inv));
        return;
    }
    loop {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let delta = (seed >> 11) & ((1u64 << k.m) - 1);
        if delta == 0 {
            continue;
        }
        // T = sum_{i<m} (delta x)^(2^i) mod g
        let mut term: Poly = trim(vec![0, delta]);
        let mut t = term.clone();
        for _ in 1..k.m {
            term = mulmod(k, &term, &term, g);
            t = add(&t, &term);
        }
        let h = gcd(k, g, &t);
        let dh = deg(&h);
        if dh > 0 && dh < d {
            let (q, _) = divrem(k, g, &h);
            split(k, &h, out, seed ^ 0x9e37);
            split(k, &q, out, seed ^ 0x7f4a);
            return;
        }
    }
}

/// x-coordinates of the rational 3-torsion kernels (roots of psi_3).
pub fn psi3_roots(k: &F2m, b: u64) -> Vec<u64> {
    roots(k, &vec![b, 0, 0, 1, 1])
}

/// The rational 3-isogenous neighbours of E_b (a2 = 1), as b-values.
pub fn neighbours(k: &F2m, b: u64) -> Vec<u64> {
    psi3_roots(k, b).into_iter().map(|x0| b ^ x0 ^ k.sqr(x0)).collect()
}

#[derive(Debug, Clone)]
pub struct Volcano3 {
    pub height: u32,
    pub root_neighbours: usize,
    /// Distance from the root to the floor (0 = on the floor).
    pub distance_to_floor: Option<u32>,
    /// Level counted from the crater (0 = crater) = height - distance.
    pub level_from_crater: Option<u32>,
    pub vertices_explored: usize,
}

/// Position of E_b in its 3-volcano of the given height (v_3 of the Frobenius
/// conductor). Floor vertices have exactly one rational 3-isogeny when the
/// height is at least 1; BFS distance to the nearest floor vertex is the
/// distance from the root to the floor.
pub fn volcano3(k: &F2m, b: u64, height: u32) -> Volcano3 {
    let n0 = neighbours(k, b);
    if height == 0 {
        return Volcano3 { height, root_neighbours: n0.len(), distance_to_floor: Some(0), level_from_crater: Some(0), vertices_explored: 1 };
    }
    let mut seen: HashMap<u64, u32> = HashMap::new();
    let mut q = VecDeque::new();
    seen.insert(b, 0);
    q.push_back((b, n0.clone()));
    let mut found = None;
    while let Some((v, nv)) = q.pop_front() {
        let dv = seen[&v];
        if nv.len() == 1 {
            found = Some(dv);
            break;
        }
        if dv >= height {
            continue;
        }
        for w in nv {
            if !seen.contains_key(&w) {
                seen.insert(w, dv + 1);
                let nw = neighbours(k, w);
                q.push_back((w, nw));
            }
        }
    }
    Volcano3 {
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
    use crate::f2m::standard_modulus;

    fn naive_order(k: &F2m, a2: u64, b: u64) -> u64 {
        // 1 (infinity) + 1 (x = 0) + sum over x != 0 of 2 * [Tr(c) = 0]
        let mut n = 2u64;
        for x in 1..(1u64 << k.m) {
            let c = x ^ a2 ^ k.mul(b, k.inv(k.sqr(x)));
            if k.trace(c) == 0 {
                n += 2;
            }
        }
        n
    }

    #[test]
    fn neighbours_are_isogenous() {
        for m in [11u32, 13] {
            let k = F2m::new(standard_modulus(m));
            let mut tested = 0;
            for b in [3u64, 5, 77, 1001, 1234, 2047, 999, 4000, 31, 600] {
                let b = b & ((1 << m) - 1);
                if b == 0 {
                    continue;
                }
                let n = naive_order(&k, 1, b);
                assert_eq!(n % 4, 2);
                for nb in neighbours(&k, b) {
                    assert_eq!(naive_order(&k, 1, nb), n, "m={m} b={b} -> {nb}");
                    assert!(neighbours(&k, nb).contains(&b), "dual isogeny missing");
                    tested += 1;
                }
            }
            assert!(tested > 0);
        }
    }
}
