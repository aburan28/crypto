//! Delfs-Galbraith (2016), "Computing isogenies between supersingular elliptic curves over F_p":
//! (1) walk in the 2-isogeny graph over F_{p^2} from each curve until the j-invariant lies in F_p;
//! (2) connect the two F_p-rational j-invariants inside the F_p-subgraph (here by Galbraith's
//! bidirectional BFS over l in {2,3}); (3) concatenate.
use super::galbraith;
use super::graph::*;
use crate::field::{Rng, Zp, Zp2};

pub struct Stats {
    pub walk_steps: usize,
    pub bfs_nodes: usize,
}

/// Non-backtracking random walk over F_{p^2} until j in F_p.
pub fn walk_to_fp(f2: &Zp2, cache: &PhiCache, j: (u64, u64), max_steps: usize, rng: &mut Rng) -> Option<Path<(u64, u64)>> {
    let mut path = Path { js: vec![j], ells: vec![] };
    let mut cur = j;
    let mut prev = None;
    for _ in 0..max_steps {
        if cur.1 == 0 {
            return Some(path);
        }
        let ns: Vec<_> = cache.get(2).neighbors(f2, cur, rng);
        let nb: Vec<_> = ns.iter().copied().filter(|&n| Some(n) != prev).collect();
        let ns = if nb.is_empty() { ns } else { nb };
        let nx = ns[rng.below(ns.len() as u64) as usize];
        prev = Some(cur);
        cur = nx;
        path.js.push(cur);
        path.ells.push(2);
    }
    if cur.1 == 0 {
        Some(path)
    } else {
        None
    }
}

pub fn delfs_galbraith(
    f2: &Zp2,
    fp: &Zp,
    cache: &PhiCache,
    j1: (u64, u64),
    j2: (u64, u64),
    max_walk: usize,
    max_bfs_nodes: usize,
    rng: &mut Rng,
) -> (Option<Path<(u64, u64)>>, Stats) {
    let mut st = Stats { walk_steps: 0, bfs_nodes: 0 };
    let Some(w1) = walk_to_fp(f2, cache, j1, max_walk, rng) else { return (None, st) };
    let Some(w2) = walk_to_fp(f2, cache, j2, max_walk, rng) else { return (None, st) };
    st.walk_steps = w1.len() + w2.len();
    let (e1, e2) = (w1.js.last().unwrap().0, w2.js.last().unwrap().0);
    let (mid, gs) = galbraith::galbraith(fp, cache, &[2, 3], e1, e2, max_bfs_nodes, rng);
    st.bfs_nodes = gs.nodes_expanded;
    let Some(mid) = mid else { return (None, st) };
    let mid2 = Path { js: mid.js.iter().map(|&a| (a, 0)).collect(), ells: mid.ells.clone() };
    (Some(w1.concat(&mid2).concat(&w2.reversed())), st)
}
