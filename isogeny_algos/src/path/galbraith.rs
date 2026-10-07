//! Galbraith (1999), "Constructing isogenies between elliptic curves over finite fields":
//! bidirectional breadth-first search in the l-isogeny graph (here on the j-line, over a set of
//! small primes l), expanding the smaller frontier a level at a time.
use super::graph::*;
use crate::field::{Field, Rng};
use std::collections::HashMap;

pub struct Stats {
    pub nodes_expanded: usize,
}

pub fn galbraith<F: Field>(
    f: &F,
    cache: &PhiCache,
    ells: &[usize],
    j1: F::E,
    j2: F::E,
    max_nodes: usize,
    rng: &mut Rng,
) -> (Option<Path<F::E>>, Stats) {
    let mut st = Stats { nodes_expanded: 0 };
    if j1 == j2 {
        return (Some(Path { js: vec![j1], ells: vec![] }), st);
    }
    let mut par: [HashMap<F::E, (F::E, usize)>; 2] = [HashMap::new(), HashMap::new()];
    let mut front: [Vec<F::E>; 2] = [vec![j1], vec![j2]];
    par[0].insert(j1, (j1, 0));
    par[1].insert(j2, (j2, 0));
    let build = |par: &[HashMap<F::E, (F::E, usize)>; 2], m: F::E| -> Path<F::E> {
        // j1 -> m
        let mut js = vec![m];
        let mut ells = vec![];
        let mut c = m;
        while c != j1 {
            let (pj, l) = par[0][&c];
            js.push(pj);
            ells.push(l);
            c = pj;
        }
        js.reverse();
        ells.reverse();
        let mut first = Path { js, ells };
        // m -> j2 (follow parents of side 2; edges are symmetric)
        let mut js2 = vec![m];
        let mut ells2 = vec![];
        let mut c = m;
        while c != j2 {
            let (pj, l) = par[1][&c];
            js2.push(pj);
            ells2.push(l);
            c = pj;
        }
        first = first.concat(&Path { js: js2, ells: ells2 });
        first
    };
    loop {
        let side = if front[0].len() <= front[1].len() { 0 } else { 1 };
        if front[side].is_empty() || par[0].len() + par[1].len() > max_nodes {
            return (None, st);
        }
        let cur = std::mem::take(&mut front[side]);
        for j in cur {
            st.nodes_expanded += 1;
            for (l, jn) in neighbors(f, cache, ells, j, rng) {
                if par[side].contains_key(&jn) {
                    continue;
                }
                par[side].insert(jn, (j, l));
                if par[1 - side].contains_key(&jn) {
                    return (Some(build(&par, jn)), st);
                }
                front[side].push(jn);
            }
        }
    }
}
