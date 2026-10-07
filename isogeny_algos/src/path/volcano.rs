//! Kohel (1996): navigation of l-isogeny volcanoes. Height from the Frobenius discriminant,
//! level/ascent/descent by comparing lengths of non-backtracking walks, and the crater walk that
//! turns this into a path-finding algorithm between two curves in the same l-volcano.
use super::graph::*;
use crate::field::{Field, Rng};

/// Height of the l-volcano of an ordinary curve with trace t over F_p (l odd): v_l(conductor of Z[pi]).
pub fn volcano_height(p: u64, t: i64, ell: u64) -> u32 {
    let d = (t as i128) * (t as i128) - 4 * (p as i128);
    let mut d = d.abs();
    let mut v = 0;
    while d % (ell as i128) == 0 && d != 0 {
        d /= ell as i128;
        v += 1;
    }
    v / 2
}

fn nbrs<F: Field>(f: &F, cache: &PhiCache, ell: usize, j: F::E, rng: &mut Rng) -> Vec<F::E> {
    cache.get(ell).neighbors(f, j, rng)
}

/// Length of the non-backtracking walk starting along the edge j -> first; Some(n) if it reaches
/// a floor vertex (exactly one neighbour) after n steps, None if not within `cap` steps.
fn walk_len<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    j: F::E,
    first: F::E,
    cap: usize,
    rng: &mut Rng,
) -> Option<usize> {
    let (mut prev, mut cur) = (j, first);
    for n in 1..=cap {
        let ns = nbrs(f, cache, ell, cur, rng);
        if ns.len() <= 1 {
            return Some(n);
        }
        let next = ns.into_iter().find(|&x| x != prev)?;
        prev = cur;
        cur = next;
    }
    None
}

/// (level, index of the unique up-neighbour if level > 0) for a vertex in a volcano of height h.
pub fn level_and_up<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    h: u32,
    j: F::E,
    rng: &mut Rng,
) -> (u32, Option<F::E>) {
    let ns = nbrs(f, cache, ell, j, rng);
    if h == 0 || ns.len() <= 1 {
        // floor (level h) vertex has a single neighbour: the up edge (if h>0)
        return (h, if h > 0 { ns.first().copied() } else { None });
    }
    let cap = 2 * h as usize + 2;
    let lens: Vec<Option<usize>> = ns
        .iter()
        .map(|&n| walk_len(f, cache, ell, j, n, cap, rng))
        .collect();
    // descending edges all reach the floor in exactly h - level steps; take the minimum finite length
    let m = lens.iter().flatten().min().copied().unwrap_or(h as usize);
    let level = h - m as u32;
    let up = if level > 0 {
        ns.iter()
            .zip(&lens)
            .find(|(_, l)| **l != Some(m))
            .map(|(n, _)| *n)
    } else {
        None
    };
    (level, up)
}

/// Ascend to the crater; returns the path (j ... crater vertex).
pub fn ascend_to_crater<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    h: u32,
    j: F::E,
    rng: &mut Rng,
) -> Path<F::E> {
    let mut js = vec![j];
    let mut cur = j;
    // A vertex of a height-h volcano needs at most h ascents; the bound turns a mis-detected
    // up-edge into a visible failure (path that ends below the crater) instead of a hang.
    for _ in 0..=h {
        let (lvl, up) = level_and_up(f, cache, ell, h, cur, rng);
        match (lvl, up) {
            (l, Some(u)) if l > 0 => {
                js.push(u);
                cur = u;
            }
            _ => break,
        }
    }
    let ells = vec![ell; js.len() - 1];
    Path { js, ells }
}

/// Descend to the floor from level `l` (first step avoids the up edge).
pub fn descend_to_floor<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    h: u32,
    j: F::E,
    rng: &mut Rng,
) -> Path<F::E> {
    let (lvl, up) = level_and_up(f, cache, ell, h, j, rng);
    let mut js = vec![j];
    if lvl == h {
        return Path { js, ells: vec![] };
    }
    let mut prev = up;
    let mut cur = j;
    loop {
        let ns = nbrs(f, cache, ell, cur, rng);
        let next = ns.into_iter().find(|&x| Some(x) != prev);
        match next {
            Some(n) => {
                js.push(n);
                prev = Some(cur);
                cur = n;
            }
            None => break,
        }
    }
    let ells = vec![ell; js.len() - 1];
    Path { js, ells }
}

/// Walk around the crater from `start` (horizontal edges only) looking for `target`.
fn crater_walk<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    h: u32,
    start: F::E,
    target: F::E,
    rng: &mut Rng,
    max: usize,
) -> Option<Path<F::E>> {
    if start == target {
        return Some(Path {
            js: vec![start],
            ells: vec![],
        });
    }
    let horizontals = |j: F::E, rng: &mut Rng| -> Vec<F::E> {
        let ns = nbrs(f, cache, ell, j, rng);
        ns.into_iter()
            .filter(|&n| {
                let (lv, _) = level_and_up(f, cache, ell, h, n, rng);
                lv == 0
            })
            .collect()
    };
    for first in horizontals(start, rng) {
        let mut js = vec![start, first];
        let (mut prev, mut cur) = (start, first);
        for _ in 0..max {
            if cur == target {
                return Some(Path {
                    ells: vec![ell; js.len() - 1],
                    js,
                });
            }
            if cur == start {
                break;
            }
            let hs = horizontals(cur, rng);
            let Some(next) = hs.into_iter().find(|&x| x != prev) else {
                break;
            };
            prev = cur;
            cur = next;
            js.push(cur);
        }
    }
    None
}

/// Kohel's algorithm: path between two curves of the same l-volcano (same crater cycle):
/// ascend both to the crater, walk the crater, descend (reversed ascent of the second).
pub fn kohel_volcano_path<F: Field>(
    f: &F,
    cache: &PhiCache,
    ell: usize,
    h: u32,
    j1: F::E,
    j2: F::E,
    rng: &mut Rng,
) -> Option<Path<F::E>> {
    let up1 = ascend_to_crater(f, cache, ell, h, j1, rng);
    let up2 = ascend_to_crater(f, cache, ell, h, j2, rng);
    let (t1, t2) = (*up1.js.last().unwrap(), *up2.js.last().unwrap());
    let crater = crater_walk(f, cache, ell, h, t1, t2, rng, 1 << 20)?;
    Some(up1.concat(&crater).concat(&up2.reversed()))
}
