//! Galbraith-Hess-Smart (2002), "Extending the GHS Weil descent attack" / fast isogeny
//! computation: random non-backtracking walks from both curves until the visited j-sets collide
//! (birthday, O(sqrt(#class)) steps). `ghs_volcano` first ascends both curves to the craters of
//! every l-volcano of positive height (Kohel), so the walks run on a single endomorphism order.
use super::graph::*;
use super::volcano;
use crate::field::{Field, Rng};
use std::collections::HashMap;

pub struct Stats {
    pub steps: usize,
}

struct Walker<E> {
    walk: Vec<E>,
    ells: Vec<usize>, // ells[i]: edge walk[i] -> walk[i+1]
    seen: HashMap<E, usize>,
}

pub fn ghs<F: Field>(
    f: &F,
    cache: &PhiCache,
    ells: &[usize],
    j1: F::E,
    j2: F::E,
    max_steps: usize,
    rng: &mut Rng,
) -> (Option<Path<F::E>>, Stats) {
    let mut st = Stats { steps: 0 };
    let mut w: [Walker<F::E>; 2] = [
        Walker {
            walk: vec![j1],
            ells: vec![],
            seen: HashMap::from([(j1, 0)]),
        },
        Walker {
            walk: vec![j2],
            ells: vec![],
            seen: HashMap::from([(j2, 0)]),
        },
    ];
    if j1 == j2 {
        return (
            Some(Path {
                js: vec![j1],
                ells: vec![],
            }),
            st,
        );
    }
    let mut stuck = [false; 2];
    let mut revisits = [0usize; 2]; // consecutive steps landing on an already-visited node
    while st.steps < max_steps && !(stuck[0] && stuck[1]) {
        for side in 0..2 {
            st.steps += 1;
            let cur = *w[side].walk.last().unwrap();
            let prev = if w[side].walk.len() >= 2 {
                Some(w[side].walk[w[side].walk.len() - 2])
            } else {
                None
            };
            let mut ns = neighbors(f, cache, ells, cur, rng);
            let non_back: Vec<_> = ns
                .iter()
                .copied()
                .filter(|&(_, n)| Some(n) != prev)
                .collect();
            if !non_back.is_empty() {
                ns = non_back;
            }
            if ns.is_empty() {
                stuck[side] = true; // no rational neighbour for any walk prime: this side cannot move
                continue;
            }
            let (l, nx) = ns[rng.below(ns.len() as u64) as usize];
            let idx = w[side].walk.len();
            w[side].walk.push(nx);
            w[side].ells.push(l);
            if w[side].seen.contains_key(&nx) {
                revisits[side] += 1;
                // saturated: every reachable node has been visited ~20 times over
                if revisits[side] > 20 * w[side].seen.len() + 100 {
                    stuck[side] = true;
                }
            } else {
                revisits[side] = 0;
                w[side].seen.insert(nx, idx);
            }
            if let Some(&k) = w[1 - side].seen.get(&nx) {
                let (a, ia, b, ib) = if side == 0 {
                    (0, idx, 1, k)
                } else {
                    (0, k, 1, idx)
                };
                let pa = Path {
                    js: w[a].walk[..=ia].to_vec(),
                    ells: w[a].ells[..ia].to_vec(),
                };
                let pb = Path {
                    js: w[b].walk[..=ib].to_vec(),
                    ells: w[b].ells[..ib].to_vec(),
                };
                return (Some(pa.concat(&pb.reversed())), st);
            }
        }
    }
    (None, st)
}

/// GHS with Kohel volcano normalisation. `trace` is the Frobenius trace of both curves.
pub fn ghs_volcano<F: Field>(
    f: &F,
    cache: &PhiCache,
    ells: &[usize],
    p: u64,
    trace: i64,
    j1: F::E,
    j2: F::E,
    max_steps: usize,
    rng: &mut Rng,
) -> (Option<Path<F::E>>, Stats) {
    let mut a = Path {
        js: vec![j1],
        ells: vec![],
    };
    let mut b = Path {
        js: vec![j2],
        ells: vec![],
    };
    let mut walk_ells = vec![];
    for &l in ells {
        let h = volcano::volcano_height(p, trace, l as u64);
        if h > 0 {
            let ua = volcano::ascend_to_crater(f, cache, l, h, *a.js.last().unwrap(), rng);
            a = a.concat(&ua);
            let ub = volcano::ascend_to_crater(f, cache, l, h, *b.js.last().unwrap(), rng);
            b = b.concat(&ub);
        } else {
            walk_ells.push(l);
        }
    }
    if walk_ells.is_empty() {
        return (None, Stats { steps: 0 });
    }
    let (mid, st) = ghs(
        f,
        cache,
        &walk_ells,
        *a.js.last().unwrap(),
        *b.js.last().unwrap(),
        max_steps,
        rng,
    );
    (mid.map(|m| a.concat(&m).concat(&b.reversed())), st)
}
