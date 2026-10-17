//! Step 1 of the fixed-degree algorithm: End(E) of an arbitrary supersingular curve E over
//! F_{p^2}, as a left O_0-ideal I with O_R(I) = End(E) (Deuring correspondence on E_0: y^2 =
//! x^3 + x, p = 3 mod 4).
//!
//! This is the classical form of the reduction of Eisenträger–Hallgren–Lauter–Morrison–Petit
//! (EUROCRYPT 2018) and Eisenträger et al. (Open Book Series 2020): computing End(E) reduces to
//! finding one isogeny path from E_0 to E, and the path is found by a generic search.
//!
//! 1. **Path.** Bidirectional breadth-first search in the 2-isogeny graph on the j-line over
//!    F_{p^2} (neighbours = distinct roots of Phi_2(j, Y)) from j = 1728 and from j(E) until the
//!    two balls meet: O~(sqrt p) vertices, the classical O~(p^{1/2}) cost quoted on the slides.
//!    (The quantum O~(p^{1/4}) variant replaces this search by Grover; it is not simulated.)
//! 2. **Translation.** Walk the path on the quaternion side: if the current left O_0-ideal I has
//!    Deuring image E_I (j(E_I) = j_t), the l-isogenies out of E_I correspond to the l + 1 left
//!    O_0-ideals I' = I l + I beta (beta in O_R(I), l | Nrd(beta), beta primitive mod l) of norm
//!    N(I) l (`sub_ideals`); the one whose Deuring image has j = j_{t+1} is the next ideal. After
//!    each step I is replaced by an equivalent ideal of minimal norm (`reduce_ideal`), which
//!    leaves its class, hence the curve, unchanged and keeps the norms near sqrt(p).
//!
//! Deuring images are computed by `represent::fast_ideal_to_j` (equivalent ideal of norm
//! dividing the odd smooth part of p^2 - 1, kernel by 2-dimensional discrete logarithms over
//! F_{p^4}, Vélu chain), so p < 2^31 and the instance must have enough odd smooth torsion.
use crate::field::{Rng, Zp2};
use crate::find::modpoly::Phi;
use crate::int::Int;
use crate::quat::deuring::Deuring;
use crate::quat::{Alg, Lattice, Quat};
use std::collections::HashMap;

pub type J = (u64, u64);

/// An equivalent left ideal of minimal norm: J = I conj(delta) / N(I) for a shortest delta in I.
/// Returns (J, xi) with J = I xi.
pub fn reduce_ideal(alg: &Alg, o0: &Lattice, i: &Lattice) -> (Lattice, Quat) {
    let (ni, nd) = i.norm(o0);
    assert!(nd == Int::one(), "reduce_ideal: integral ideals only");
    let rb = i.reduced_basis(alg);
    let v = rb.iter().min_by(|a, b| alg.qf(a).cmp(&alg.qf(b))).unwrap();
    let delta = Quat::new(v.clone(), i.den.clone());
    let xi = delta.conj().scale(&Int::one(), &ni);
    (i.rmul(alg, &xi), xi)
}

/// The l + 1 left O_0-ideals I' of norm N(I) l contained in I (I' = I l + I beta, beta in
/// O_R(I) primitive mod l with l | Nrd(beta)); these are the ideals of the l-isogenies out of
/// the Deuring image of I. Enumerates beta over O_R(I) / l O_R(I) (l^4 classes).
pub fn sub_ideals(alg: &Alg, o0: &Lattice, i: &Lattice, ell: u64) -> Vec<Lattice> {
    let ni = i.norm_int(o0);
    let target = &ni * &Int::from(ell);
    let or = i.right_order(alg, o0);
    let ob = or.basis();
    let il: Vec<Quat> = i
        .basis()
        .iter()
        .map(|b| b.scale(&Int::from(ell), &Int::one()))
        .collect();
    let mut out: Vec<Lattice> = vec![];
    let l = ell as usize;
    for code in 1..l.pow(4) {
        let mut c = [0i64; 4];
        let mut t = code;
        for slot in c.iter_mut() {
            *slot = (t % l) as i64;
            t /= l;
        }
        let mut beta = Quat::zero();
        for (k, b) in ob.iter().enumerate() {
            if c[k] != 0 {
                beta = beta.add(&b.scale(&Int::from(c[k]), &Int::one()));
            }
        }
        let (n, d) = alg.nrd(&beta);
        if d != Int::one() || !n.modulo(&Int::from(ell)).is_zero() {
            continue;
        }
        let mut gens = il.clone();
        gens.extend(i.rmul(alg, &beta).basis());
        let cand = Lattice::from_gens(&gens);
        if cand.norm(o0) != (target.clone(), Int::one()) {
            continue;
        }
        if !out.contains(&cand) {
            out.push(cand);
            if out.len() == l + 1 {
                break;
            }
        }
    }
    out
}

/// Distinct neighbours of j in the l-isogeny graph on the j-line (roots of Phi_l(j, Y)).
pub struct JGraph {
    pub f2: Zp2,
    pub phi: Phi<Zp2>,
    pub cache: HashMap<J, Vec<J>>,
    rng: Rng,
}

impl JGraph {
    pub fn new(p: u64, ell: usize) -> JGraph {
        let fp = crate::field::Zp::new(p);
        JGraph {
            f2: Zp2::new(p),
            phi: Phi::compute(&fp, ell).lift(|x| (x, 0u64)),
            cache: HashMap::new(),
            rng: Rng::new(0x9a7f),
        }
    }
    pub fn neighbours(&mut self, j: J) -> Vec<J> {
        if let Some(v) = self.cache.get(&j) {
            return v.clone();
        }
        let mut v = self.phi.neighbors(&self.f2, j, &mut self.rng);
        v.sort();
        v.dedup();
        self.cache.insert(j, v.clone());
        v
    }
}

/// Bidirectional BFS: a shortest j-path from `a` to `b` (inclusive), and the number of
/// vertices visited.
pub fn bfs_path(g: &mut JGraph, a: J, b: J) -> Option<(Vec<J>, usize)> {
    if a == b {
        return Some((vec![a], 1));
    }
    let mut par: [HashMap<J, Option<J>>; 2] = [HashMap::new(), HashMap::new()];
    par[0].insert(a, None);
    par[1].insert(b, None);
    let mut front: [Vec<J>; 2] = [vec![a], vec![b]];
    loop {
        if front[0].is_empty() && front[1].is_empty() {
            return None;
        }
        // expand the smaller nonempty frontier by one layer
        let s = if front[1].is_empty() || (!front[0].is_empty() && front[0].len() <= front[1].len())
        {
            0
        } else {
            1
        };
        let mut next = vec![];
        let mut meet = None;
        for &u in &front[s] {
            for v in g.neighbours(u) {
                if par[s].contains_key(&v) {
                    continue;
                }
                par[s].insert(v, Some(u));
                if par[1 - s].contains_key(&v) && meet.is_none() {
                    meet = Some(v);
                }
                next.push(v);
            }
        }
        front[s] = next;
        if let Some(m) = meet {
            let walk = |side: usize, par: &[HashMap<J, Option<J>>; 2]| {
                let mut out = vec![m];
                let mut cur = m;
                while let Some(Some(pv)) = par[side].get(&cur) {
                    out.push(*pv);
                    cur = *pv;
                }
                out
            };
            let mut left = walk(0, &par);
            left.reverse();
            let right = walk(1, &par);
            left.extend_from_slice(&right[1..]);
            let visited = par[0].len() + par[1].len();
            return Some((left, visited));
        }
    }
}

#[derive(Clone, Debug)]
pub struct EndRing {
    /// left O_0-ideal with O_R(ideal) = End(E) (reduced: minimal norm in its class)
    pub ideal: Lattice,
    /// j-path from 1728 to j(E) used for the translation
    pub path: Vec<J>,
    /// vertices visited by the bidirectional search
    pub graph_vertices: usize,
    /// Deuring evaluations (ideal -> j) spent in the translation
    pub deuring_evals: usize,
}

/// Translate a j-path starting at 1728 into the left O_0-ideal of its endpoint.
pub fn translate_path(dg: &Deuring, path: &[J], ell: u64) -> Option<(Lattice, usize)> {
    let (alg, o0) = (&dg.alg, &dg.o0);
    let mut ideal = o0.clone();
    let mut evals = 0usize;
    for next in &path[1..] {
        let cands = sub_ideals(alg, o0, &ideal, ell);
        let mut found = None;
        for c in cands {
            evals += 1;
            if super::represent::fast_ideal_to_j(dg, &c)? == *next {
                found = Some(c);
                break;
            }
        }
        ideal = reduce_ideal(alg, o0, &found?).0;
    }
    Some((ideal, evals))
}

/// End(E) for the supersingular curve with j-invariant `j` (Step 1 of the slides).
pub fn end_ring(dg: &Deuring, graph: &mut JGraph, j: J) -> Option<EndRing> {
    let j0 = (1728 % dg.p, 0u64);
    let (path, graph_vertices) = bfs_path(graph, j0, j)?;
    let ell = graph.phi.ell as u64;
    let (ideal, deuring_evals) = translate_path(dg, &path, ell)?;
    Some(EndRing {
        ideal,
        path,
        graph_vertices,
        deuring_evals,
    })
}

/// A random supersingular j-invariant: endpoint of a random non-backtracking 2-walk of the given
/// length from 1728 (test-instance generator; independent of the quaternion side).
pub fn random_ss_j(g: &mut JGraph, p: u64, steps: usize, rng: &mut Rng) -> J {
    let mut prev = None;
    let mut cur = (1728 % p, 0u64);
    for _ in 0..steps {
        let nb: Vec<J> = g
            .neighbours(cur)
            .into_iter()
            .filter(|&v| Some(v) != prev)
            .collect();
        if nb.is_empty() {
            break;
        }
        let nxt = nb[rng.below(nb.len() as u64) as usize];
        prev = Some(cur);
        cur = nxt;
    }
    cur
}
