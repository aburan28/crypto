//! Couveignes (2006), "Hard homogeneous spaces" / Rostovtsev-Stolbunov: ideal classes of the
//! maximal order act on the ordinary curves with that endomorphism ring. Given E1, E2 = [a]E1,
//! recover an exponent vector over a set of split primes by meet-in-the-middle over boxes.
//!
//! The direction of each prime ideal is fixed canonically by the Frobenius eigenvalue (up to sign,
//! hence twist-invariant) on the kernel of the horizontal isogeny; computing it needs the kernel
//! polynomial (from Elkies/BMSS) and arithmetic in F_p[x]/h.
use super::graph::*;
use crate::curve::*;
use crate::field::{Rng, Zp};
use crate::find::{divpoly::division_poly, elkies};
use crate::poly::{self, Poly};
use std::collections::HashMap;

/// min(lambda, l - lambda) for the Frobenius eigenvalue on the kernel with polynomial h.
pub fn eigen_class(fp: &Zp, e: &Curve<u64>, h: &Poly<Zp>, ell: u64) -> Option<u64> {
    let x = poly::x_poly(fp);
    let xp = poly::powmod(fp, &x, fp.p as u128, h);
    if xp == poly::rem(fp, &x, h) {
        return Some(1);
    }
    let fx: Poly<Zp> = vec![e.b, e.a, 0, 1];
    let fxm = poly::rem(fp, &fx, h);
    let g = |n: usize| poly::rem(fp, &division_poly(fp, e, n), h);
    for lam in 2..=((ell - 1) / 2) as usize {
        let (gl, gm, gp) = (g(lam), g(lam - 1), g(lam + 1));
        let inv = poly::invmod(fp, &poly::mulmod(fp, &gl, &gl, h), h)?;
        let num = poly::mulmod(fp, &gm, &gp, h);
        let mut ratio = poly::mulmod(fp, &num, &inv, h);
        let four_f = poly::scale(fp, &fxm, 4);
        if lam % 2 == 0 {
            // x_lam = x - g_{l-1} g_{l+1} / (4 F g_l^2)
            let i4f = poly::invmod(fp, &four_f, h)?;
            ratio = poly::mulmod(fp, &ratio, &i4f, h);
        } else {
            ratio = poly::mulmod(fp, &ratio, &four_f, h);
        }
        let xl = poly::sub(fp, &poly::rem(fp, &x, h), &ratio);
        if xl == xp {
            return Some(lam as u64);
        }
    }
    None
}

pub struct Action<'a> {
    pub fp: &'a Zp,
    pub cache: &'a PhiCache<Zp>,
    pub plus_class: HashMap<usize, u64>,
}

impl<'a> Action<'a> {
    /// The two horizontal neighbours of j for prime l with their eigenvalue classes.
    fn oriented(&self, j: u64, ell: usize, rng: &mut Rng) -> Option<Vec<(u64, u64)>> {
        if j == 0 || j == 1728 {
            return None;
        }
        let e = from_j(self.fp, j);
        let isos = elkies::elkies_isogenies(self.fp, self.cache.get(ell), &e, rng);
        if isos.len() != 2 {
            return None;
        }
        let mut out = vec![];
        for iso in isos {
            let c = eigen_class(self.fp, &e, &iso.ker, ell as u64)?;
            out.push((jinv(self.fp, &iso.cod), c));
        }
        if out[0].1 == out[1].1 {
            return None;
        }
        Some(out)
    }
    /// Fix the + orientation of prime l at the base curve (first neighbour is +).
    pub fn init_prime(&mut self, j: u64, ell: usize, rng: &mut Rng) -> bool {
        match self.oriented(j, ell, rng) {
            Some(o) => {
                self.plus_class.insert(ell, o[0].1);
                true
            }
            None => false,
        }
    }
    pub fn step(&self, j: u64, ell: usize, positive: bool, rng: &mut Rng) -> Option<u64> {
        let o = self.oriented(j, ell, rng)?;
        let pc = self.plus_class[&ell];
        o.iter()
            .find(|&&(_, c)| (c == pc) == positive)
            .map(|&(jj, _)| jj)
    }
    pub fn apply(
        &self,
        j: u64,
        primes: &[usize],
        exps: &[i64],
        rng: &mut Rng,
    ) -> Option<Path<u64>> {
        let mut path = Path {
            js: vec![j],
            ells: vec![],
        };
        let mut cur = j;
        for (i, &l) in primes.iter().enumerate() {
            for _ in 0..exps[i].abs() {
                cur = self.step(cur, l, exps[i] > 0, rng)?;
                path.js.push(cur);
                path.ells.push(l);
            }
        }
        Some(path)
    }
}

pub struct Stats {
    pub nodes: usize,
}

/// Meet-in-the-middle over boxes [0,m]^k from both curves. Returns the exponent vector e with
/// [e]E1 = E2 and the j-path.
pub fn couveignes<'a>(
    act: &Action<'a>,
    primes: &[usize],
    m: usize,
    j1: u64,
    j2: u64,
    rng: &mut Rng,
) -> (Option<(Vec<i64>, Path<u64>)>, Stats) {
    let k = primes.len();
    let base = m + 1;
    let total = base.pow(k as u32);
    let mut st = Stats { nodes: 0 };
    let grow = |start: u64,
                rng: &mut Rng,
                st: &mut Stats,
                other: Option<&HashMap<u64, usize>>|
     -> (Vec<u64>, Option<(usize, usize)>) {
        let mut nodes = vec![start; total];
        let mut hit = None;
        for n in 0..total {
            if n > 0 {
                let mut digit = n;
                let mut i = 0;
                while digit % base == 0 {
                    digit /= base;
                    i += 1;
                }
                let parent = n - base.pow(i as u32);
                if nodes[parent] == u64::MAX {
                    nodes[n] = u64::MAX;
                    continue;
                }
                match act.step(nodes[parent], primes[i], true, rng) {
                    Some(j) => nodes[n] = j,
                    None => {
                        nodes[n] = u64::MAX;
                        continue;
                    }
                }
                st.nodes += 1;
            }
            if let Some(o) = other {
                if let Some(&a) = o.get(&nodes[n]) {
                    hit = Some((a, n));
                    break;
                }
            }
        }
        (nodes, hit)
    };
    let (nodes_a, _) = grow(j1, rng, &mut st, None);
    let mut map_a: HashMap<u64, usize> = HashMap::new();
    for (i, &j) in nodes_a.iter().enumerate() {
        if j != u64::MAX {
            map_a.entry(j).or_insert(i);
        }
    }
    let (_nodes_b, hit) = grow(j2, rng, &mut st, Some(&map_a));
    let Some((na, nb)) = hit else {
        return (None, st);
    };
    let digits = |mut n: usize| -> Vec<i64> {
        (0..k)
            .map(|_| {
                let d = (n % base) as i64;
                n /= base;
                d
            })
            .collect()
    };
    let (a, b) = (digits(na), digits(nb));
    let e: Vec<i64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
    let pa = act.apply(j1, primes, &a, rng).unwrap();
    let pb = act.apply(j2, primes, &b, rng).unwrap();
    (Some((e, pa.concat(&pb.reversed()))), st)
}

/// Primes (from the cached set) that split with distinguishable orientation at j.
pub fn select_primes(act: &mut Action, j: u64, rng: &mut Rng) -> Vec<usize> {
    let mut out = vec![];
    for l in act.cache.ells() {
        if l >= 3 && act.init_prime(j, l, rng) {
            out.push(l);
        }
    }
    out
}
