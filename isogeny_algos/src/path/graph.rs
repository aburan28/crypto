//! The l-isogeny graph on the j-line, plus conversion of a j-path into explicit isogenies.
use crate::curve::*;
use crate::field::{Field, Rng};
use crate::find::{elkies, modpoly::Phi};
use std::collections::HashMap;

pub struct PhiCache<F: Field> {
    pub phis: Vec<Phi<F>>,
}
impl<F: Field> PhiCache<F> {
    pub fn new(f: &F, ells: &[usize]) -> Self {
        PhiCache {
            phis: ells.iter().map(|&l| Phi::compute(f, l)).collect(),
        }
    }
    /// The cache over another field through an embedding (e.g. F_p into F_{p^2}).
    pub fn lift<G: Field>(&self, emb: impl Fn(F::E) -> G::E + Copy) -> PhiCache<G> {
        PhiCache {
            phis: self.phis.iter().map(|p| p.lift(emb)).collect(),
        }
    }
    pub fn get(&self, ell: usize) -> &Phi<F> {
        self.phis
            .iter()
            .find(|p| p.ell == ell)
            .expect("Phi not cached")
    }
    pub fn ells(&self) -> Vec<usize> {
        self.phis.iter().map(|p| p.ell).collect()
    }
}

#[derive(Clone, Debug)]
pub struct Path<E> {
    pub js: Vec<E>,
    pub ells: Vec<usize>,
}
impl<E: Copy> Path<E> {
    pub fn len(&self) -> usize {
        self.ells.len()
    }
    pub fn reversed(&self) -> Path<E> {
        Path {
            js: self.js.iter().rev().copied().collect(),
            ells: self.ells.iter().rev().copied().collect(),
        }
    }
    pub fn concat(&self, other: &Path<E>) -> Path<E> {
        let mut js = self.js.clone();
        js.extend_from_slice(&other.js[1..]);
        let mut ells = self.ells.clone();
        ells.extend_from_slice(&other.ells);
        Path { js, ells }
    }
}

/// A neighbour oracle on some isogeny graph: the (l, j') edges out of j. Path-finding algorithms
/// (Galbraith BFS, GHS walks) are written against this, so the same code runs with modular
/// polynomials over F_p, F_{p^2} or GF(2^n), or with kernel-polynomial factoring.
pub trait Oracle<E> {
    fn neighbors(&self, j: E, rng: &mut Rng) -> Vec<(usize, E)>;
}

/// Roots of Phi_l(j, Y) for each l in `ells`.
pub struct PhiOracle<'a, F: Field> {
    pub f: &'a F,
    pub cache: &'a PhiCache<F>,
    pub ells: &'a [usize],
}

impl<F: Field> Oracle<F::E> for PhiOracle<'_, F> {
    fn neighbors(&self, j: F::E, rng: &mut Rng) -> Vec<(usize, F::E)> {
        neighbors(self.f, self.cache, self.ells, j, rng)
    }
}

/// Distinct F-rational neighbours of j in the l-graph for each l in `ells`.
pub fn neighbors<F: Field>(
    f: &F,
    cache: &PhiCache<F>,
    ells: &[usize],
    j: F::E,
    rng: &mut Rng,
) -> Vec<(usize, F::E)> {
    let mut out = vec![];
    for &l in ells {
        for jn in cache.get(l).neighbors(f, j, rng) {
            out.push((l, jn));
        }
    }
    out
}

pub fn verify_path<F: Field>(f: &F, cache: &PhiCache<F>, path: &Path<F::E>) -> bool {
    path.ells
        .iter()
        .enumerate()
        .all(|(i, &l)| cache.get(l).eval(f, path.js[i], path.js[i + 1]) == f.zero())
}

/// Turn a j-path into explicit normalised isogenies starting at `e`.
/// Fails (None) if some step passes through j in {0,1728} or a degenerate Phi derivative.
pub fn explicit_chain<F: Field>(
    f: &F,
    cache: &PhiCache<F>,
    e: &Curve<F::E>,
    path: &Path<F::E>,
) -> Option<Vec<RatIsogeny<F>>> {
    let mut cur = *e;
    if jinv(f, &cur) != path.js[0] {
        return None;
    }
    let mut out = vec![];
    for (i, &l) in path.ells.iter().enumerate() {
        let phi = cache.get(l);
        let ep = elkies::elkies_codomain(f, phi, &cur, path.js[i + 1])?;
        let iso = elkies::bmss_isogeny(f, &cur, &ep, l)?;
        cur = iso.cod;
        out.push(iso);
    }
    Some(out)
}

/// Apply a chain to a point.
pub fn apply_chain<F: Field>(f: &F, chain: &[RatIsogeny<F>], p: &Pt<F::E>) -> Pt<F::E> {
    let mut q = *p;
    for iso in chain {
        q = iso.eval(f, &q);
    }
    q
}

pub fn index_map<E: std::hash::Hash + Eq + Copy>(v: &[E]) -> HashMap<E, usize> {
    let mut m = HashMap::new();
    for (i, &x) in v.iter().enumerate() {
        m.entry(x).or_insert(i);
    }
    m
}
