//! Baby-step giant-step search for every N in a Hasse window with N P = O.

use std::collections::HashMap;

pub trait Group {
    type P: Copy + Eq + std::fmt::Debug;
    fn inf(&self) -> Self::P;
    fn is_inf(&self, p: Self::P) -> bool;
    fn add(&self, p: Self::P, q: Self::P) -> Self::P;
    fn neg(&self, p: Self::P) -> Self::P;
    fn mul(&self, p: Self::P, n: u128) -> Self::P;
    /// A key shared by exactly {P, -P} (the x-coordinate).
    fn key(&self, p: Self::P) -> u64;
}

/// All N in [lo, hi] with N = c (mod k) and N P = O. `None` when the baby
/// table meets O (the stride k P has order <= the table size).
pub fn multiples_in_window<G: Group>(
    g: &G,
    p: G::P,
    lo: u128,
    hi: u128,
    k: u128,
    c: u128,
) -> Option<Vec<u128>> {
    let n0 = lo + ((c + k - lo % k) % k);
    if n0 > hi {
        return Some(vec![]);
    }
    let span = (hi - n0) / k; // i in [0, span]
    let r = g.mul(p, k);
    if g.is_inf(r) {
        return None;
    }
    let q = g.mul(p, n0);
    let s = (((span / 2) as f64).sqrt() as u128 + 1).max(1);
    let mut table: HashMap<u64, (u32, G::P)> = HashMap::with_capacity(s as usize * 2);
    let mut jr = r;
    for j in 1..=s {
        if g.is_inf(jr) {
            return None;
        }
        table.entry(g.key(jr)).or_insert((j as u32, jr));
        jr = g.add(jr, r);
    }
    let step = 2 * s;
    let giant = g.mul(r, step);
    let mut t = q;
    let mut out = Vec::new();
    let mut base: u128 = 0;
    loop {
        if base > span + s {
            break;
        }
        if g.is_inf(t) {
            if base <= span {
                out.push(base);
            }
        } else if let Some(&(j, pt)) = table.get(&g.key(t)) {
            let j = j as u128;
            if pt == t {
                // Q + base R = j R  =>  i = base - j
                if base >= j && base - j <= span {
                    out.push(base - j);
                }
            } else {
                // Q + base R = -j R  =>  i = base + j
                debug_assert_eq!(g.neg(pt), t);
                if base + j <= span {
                    out.push(base + j);
                }
            }
        }
        t = g.add(t, giant);
        base += step;
    }
    let mut res: Vec<u128> = out.into_iter().map(|i| n0 + k * i).collect();
    res.sort_unstable();
    res.dedup();
    // A key collision between distinct baby multiples can hide a match only if
    // two baby entries share an x-coordinate, i.e. j R = -j' R with j != j',
    // which means ord(R) <= 2s; reject that case explicitly.
    if table.len() < s as usize {
        return None;
    }
    for &n in &res {
        debug_assert!(g.is_inf(g.mul(p, n)));
    }
    Some(res)
}
