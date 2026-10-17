//! Volcano position and 2-sum factor-base yields of ordinary binary curves
//! y^2 + xy = x^3 + a2 x^2 + a6 over GF(2^n) (`binary.rs`): the Frobenius discriminant split
//! D = f_pi^2 D_K, ell-isogeny neighbours with multiplicity (roots of Phi_ell mod 2), Kohel
//! levels by non-backtracking walks, the connected component of a curve in the ell-graph with a
//! vertex bound, and the exact 2-sum relation yield of the factor base
//! V = span_F2{1, z, ..., z^(l-1)} = {x < 2^l} in the polynomial basis.
use crate::binary::{BPt, BinCurve};
use crate::curve::Pt;
use crate::field::{factor_u64, is_prime, Field, Rng};
use crate::find::modpoly::Phi;
use crate::gf2n::GF2n;
use crate::poly::{self, Poly};
use std::collections::HashMap;

/// D = t^2 - 4q = f_pi^2 D_K for a curve of the given order over F_q.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobDisc {
    pub trace: i128,
    pub d: i128,
    /// fundamental discriminant of Q(pi)
    pub d_k: i128,
    /// conductor of Z[pi] in the maximal order
    pub f_pi: u128,
    /// prime factorisation of f_pi (complete when `unfactored` is None)
    pub f_pi_factors: Vec<(u64, u32)>,
    /// part of |D| that was not factored (only possible above 2^64), if any
    pub unfactored: Option<u128>,
}

fn isqrt_u128(n: u128) -> u128 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u128;
    while x * x > n {
        x -= 1;
    }
    while (x + 1) * (x + 1) <= n {
        x += 1;
    }
    x
}

/// Prime factorisation of n: complete (trial division + Pollard rho) below 2^64; above, trial
/// division to 10^6, a perfect-square test and a Pollard rho on the remainder when it fits a
/// u64, otherwise the remainder is returned unfactored.
pub fn factor_abs(n: u128) -> (Vec<(u64, u32)>, Option<u128>) {
    if n <= u64::MAX as u128 {
        return (factor_u64(n as u64), None);
    }
    let mut m = n;
    let mut out: Vec<(u64, u32)> = vec![];
    let mut p = 2u64;
    while p <= 1_000_000 && (p as u128) * (p as u128) <= m {
        if m % (p as u128) == 0 {
            let mut e = 0;
            while m % (p as u128) == 0 {
                m /= p as u128;
                e += 1;
            }
            out.push((p, e));
        }
        p += if p == 2 { 1 } else { 2 };
    }
    let mut unfactored = None;
    if m > 1 {
        let s = isqrt_u128(m);
        if s * s == m && s <= u64::MAX as u128 && is_prime(s as u64) {
            out.push((s as u64, 2));
        } else if m <= u64::MAX as u128 {
            for (q, e) in factor_u64(m as u64) {
                out.push((q, e));
            }
        } else {
            unfactored = Some(m);
        }
    }
    out.sort();
    (out, unfactored)
}

/// Frobenius discriminant of a curve with `order` points over F_q and its decomposition.
pub fn frobenius_disc(q: u128, order: u128) -> FrobDisc {
    let trace = q as i128 + 1 - order as i128;
    let d = trace * trace - 4 * q as i128;
    assert!(d < 0, "not an elliptic curve order in the Hasse interval");
    let (fac, unfactored) = factor_abs(d.unsigned_abs());
    let mut f: u128 = 1;
    let mut ff: Vec<(u64, u32)> = vec![];
    for &(p, e) in &fac {
        if e >= 2 {
            f *= (p as u128).pow(e / 2);
            ff.push((p, e / 2));
        }
    }
    let mut d_k = d / (f * f) as i128;
    if d_k.rem_euclid(4) > 1 {
        // d_k = 2, 3 mod 4: D_K = 4 d_k and f loses one factor 2 (D = 0, 1 mod 4 forces 2 | f)
        d_k *= 4;
        f /= 2;
        if let Some(pos) = ff.iter().position(|&(p, _)| p == 2) {
            ff[pos].1 -= 1;
            if ff[pos].1 == 0 {
                ff.remove(pos);
            }
        }
    }
    FrobDisc {
        trace,
        d,
        d_k,
        f_pi: f,
        f_pi_factors: ff,
        unfactored,
    }
}

impl FrobDisc {
    /// v_ell(f_pi): the height of the ell-volcano.
    pub fn height(&self, ell: u64) -> u32 {
        self.f_pi_factors
            .iter()
            .find(|&&(p, _)| p == ell)
            .map_or(0, |&(_, e)| e)
    }
    /// Odd primes ell <= ell_max dividing f_pi, with their heights.
    pub fn odd_heights(&self, ell_max: u64) -> Vec<(u64, u32)> {
        self.f_pi_factors
            .iter()
            .filter(|&&(p, _)| p % 2 == 1 && p <= ell_max)
            .map(|&(p, e)| (p, e))
            .collect()
    }
}

/// Kronecker symbol (d / ell) for an odd prime ell.
pub fn kronecker_odd(d: i128, ell: u64) -> i32 {
    let a = d.rem_euclid(ell as i128) as u64;
    if a == 0 {
        return 0;
    }
    let mut r = 1u64;
    let (mut b, mut e) = (a % ell, (ell - 1) / 2);
    while e > 0 {
        if e & 1 == 1 {
            r = (r as u128 * b as u128 % ell as u128) as u64;
        }
        b = (b as u128 * b as u128 % ell as u128) as u64;
        e >>= 1;
    }
    if r == 1 {
        1
    } else {
        -1
    }
}

/// Odd primes ell <= ell_max for which 4q is a square modulo ell, i.e. the only odd ell that can
/// divide t^2 - 4q at all.
pub fn admissible_ells(q: u128, ell_max: u64) -> Vec<u64> {
    (3..=ell_max)
        .filter(|&l| l % 2 == 1 && is_prime(l))
        .filter(|&l| kronecker_odd((4 * q) as i128, l) >= 0)
        .collect()
}

/// An element of absolute trace 1 (1 itself when n is odd): the curves with a2 = 0 and a2 = this
/// element are the two quadratic twists of each j.
pub fn trace_one_element(f: &GF2n) -> u64 {
    (0..f.n)
        .map(|i| 1u64 << i)
        .find(|&b| f.trace(b) == 1)
        .expect("some basis element has trace 1")
}

/// Roots of `p` in GF(2^n) with multiplicity (repeated division by Y + r), sorted.
pub fn roots_with_multiplicity(f: &GF2n, p: &Poly<GF2n>, rng: &mut Rng) -> Vec<u64> {
    let mut out = vec![];
    let mut rs = poly::roots(f, p, rng);
    rs.sort_unstable();
    let mut cur = poly::monic(f, p);
    for r in rs {
        let lin = vec![r, 1u64]; // Y + r
        loop {
            let (q, rem) = poly::divrem(f, &cur, &lin);
            if !rem.is_empty() {
                break;
            }
            out.push(r);
            cur = q;
            if cur.len() <= 1 {
                break;
            }
        }
    }
    out.sort_unstable();
    out
}

/// The ell-isogeny graph on ordinary j-invariants over GF(2^n) (edges with multiplicity: the
/// roots of Phi_ell(j, Y) in GF(2^n)), restricted to one volcano of the given height.
pub struct EllGraph<'a> {
    pub f: &'a GF2n,
    pub phi: Phi<GF2n>,
    pub ell: u64,
    pub height: u32,
}

impl<'a> EllGraph<'a> {
    pub fn new(f: &'a GF2n, ell: u64, height: u32) -> Self {
        EllGraph {
            f,
            phi: crate::binary::phi_mod2(f, ell as usize),
            ell,
            height,
        }
    }
    /// Rational ell-isogenies out of j, as target j-invariants with multiplicity, sorted.
    pub fn neighbors(&self, j: u64, rng: &mut Rng) -> Vec<u64> {
        roots_with_multiplicity(self.f, &self.phi.y_poly(self.f, j), rng)
    }
    /// Length of the non-backtracking walk j -> first -> ... if it reaches a floor vertex
    /// (one rational ell-isogeny) within `height` steps.
    fn walk_to_floor(&self, j: u64, first: u64, rng: &mut Rng) -> Option<usize> {
        let (mut prev, mut cur) = (j, first);
        for n in 1..=self.height as usize {
            let mut ns = self.neighbors(cur, rng);
            if ns.len() <= 1 {
                return Some(n);
            }
            // drop the edge we came along (one copy: double edges and loops keep the other)
            if let Some(pos) = ns.iter().position(|&x| x == prev) {
                ns.remove(pos);
            }
            let next = *ns.first()?;
            prev = cur;
            cur = next;
        }
        None
    }
    /// Kohel level of j (0 = crater, height = floor): height minus the length of the shortest
    /// non-backtracking walk to a floor vertex, which is the straight descent. None when no walk
    /// of at most `height` steps reaches the floor (j not in a volcano of this height).
    pub fn level(&self, j: u64, rng: &mut Rng) -> Option<u32> {
        let h = self.height;
        if h == 0 {
            return Some(0);
        }
        let ns = self.neighbors(j, rng);
        if ns.len() <= 1 {
            return Some(h);
        }
        let m = ns
            .iter()
            .filter_map(|&n| self.walk_to_floor(j, n, rng))
            .min()?;
        Some(h - m as u32)
    }
}

/// Which of the two twists (a2 = pair.0, a2 = pair.1, a4 = 0, a6 = 1/j) has the given order:
/// random points are tested against the order and the twist order (the twists' orders differ
/// for an ordinary curve). None if a point is annihilated by neither (j not in the class) or the
/// tests stay ambiguous.
pub fn same_order_a2(
    f: &GF2n,
    j: u64,
    order: u128,
    a2_pair: (u64, u64),
    rng: &mut Rng,
) -> Option<u64> {
    let twist_order = 2 * f.size() + 2 - order;
    let e0 = BinCurve::new(a2_pair.0, 0, f.inv(j));
    for _ in 0..32 {
        let p = e0.random_point(f, rng);
        let z0 = e0.mul(f, &p, order) == Pt::Inf;
        let z1 = e0.mul(f, &p, twist_order) == Pt::Inf;
        match (z0, z1) {
            (true, false) => return Some(a2_pair.0),
            (false, true) => return Some(a2_pair.1),
            (false, false) => return None,
            (true, true) => {}
        }
    }
    None
}

#[derive(Clone, Debug)]
pub struct Vertex {
    pub j: u64,
    /// a2 of the model with the class order (a4 = 0, a6 = 1/j); None when `same_order_a2` failed
    pub a2: Option<u64>,
    pub parent: Option<usize>,
    pub level: Option<u32>,
    /// rational ell-isogenies out of j (neighbours with multiplicity)
    pub rational: usize,
    /// distinct neighbours
    pub neighbors: Vec<u64>,
}

#[derive(Clone, Debug)]
pub struct Component {
    pub vertices: Vec<Vertex>,
    /// the vertex bound stopped the enumeration
    pub truncated: bool,
}

/// Breadth-first enumeration of the connected component of j0 in the ell-graph, at most
/// `max_vertices` vertices, with the level and the same-order twist of every vertex.
pub fn enumerate_component(
    g: &EllGraph,
    j0: u64,
    a2_0: u64,
    order: u128,
    a2_pair: (u64, u64),
    max_vertices: usize,
    rng: &mut Rng,
    mut progress: impl FnMut(usize, usize),
) -> Component {
    let mut index: HashMap<u64, usize> = HashMap::new();
    let mut vertices: Vec<Vertex> = vec![];
    let mut truncated = false;
    index.insert(j0, 0);
    vertices.push(Vertex {
        j: j0,
        a2: Some(a2_0),
        parent: None,
        level: None,
        rational: 0,
        neighbors: vec![],
    });
    let mut head = 0;
    while head < vertices.len() {
        let j = vertices[head].j;
        let ns = g.neighbors(j, rng);
        let mut distinct = ns.clone();
        distinct.dedup();
        vertices[head].rational = ns.len();
        vertices[head].neighbors = distinct.clone();
        vertices[head].level = g.level(j, rng);
        for &jn in &distinct {
            if index.contains_key(&jn) {
                continue;
            }
            if vertices.len() >= max_vertices {
                truncated = true;
                continue;
            }
            let a2 = same_order_a2(g.f, jn, order, a2_pair, rng);
            index.insert(jn, vertices.len());
            vertices.push(Vertex {
                j: jn,
                a2,
                parent: Some(head),
                level: None,
                rational: 0,
                neighbors: vec![],
            });
        }
        head += 1;
        progress(head, vertices.len());
    }
    Component {
        vertices,
        truncated,
    }
}

/// The factor base of a curve with a4 = 0: one point (x, y) per abscissa x < 2^l on the curve
/// (x = 0 is the 2-torsion point, with its single y; every other abscissa stands for +-P).
#[derive(Clone, Debug)]
pub struct FactorBase {
    pub l: u32,
    pub xs: Vec<u64>,
    pub ys: Vec<u64>,
}

pub fn factor_base(f: &GF2n, e: &BinCurve, l: u32) -> FactorBase {
    assert!(e.a4 == 0, "normalise first");
    assert!(l <= f.n);
    let (mut xs, mut ys) = (vec![], vec![]);
    for x in 0..(1u64 << l) {
        if x == 0 {
            xs.push(0);
            ys.push(f.sqrt(e.a6).unwrap());
            continue;
        }
        let c = f.mul(e.rhs(f, x), f.inv(f.sq(x)));
        if let Some(z) = f.solve_quadratic(c) {
            xs.push(x);
            ys.push(f.mul(x, z));
        }
    }
    FactorBase { l, xs, ys }
}

impl FactorBase {
    pub fn abscissae(&self) -> usize {
        self.xs.len()
    }
    /// |F|: points with both y, the 2-torsion point once.
    pub fn points(&self) -> usize {
        2 * self.xs.len() - usize::from(self.xs.first() == Some(&0))
    }
    pub fn point(&self, i: usize) -> BPt {
        Pt::Aff(self.xs[i], self.ys[i])
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TwoSumYield {
    /// unordered pairs of distinct abscissae {x1, x2} of the factor base (P1 != +-P2)
    pub pairs: u64,
    /// distinct x(P1 + P2), x(P1 - P2) over those pairs
    pub count: u64,
    /// the same count over the pairs whose sum (resp. difference) lies in the order-r subgroup,
    /// i.e. [r]P1 = -[r]P2 (resp. [r]P1 = [r]P2); None without r
    pub subgroup_count: Option<u64>,
    /// number of (unordered pair, sign) combinations with [r]P1 + s[r]P2 = O
    pub subgroup_combos: Option<u64>,
    /// ordered (F1, F2, s), F1, F2 in F without the 2-torsion point, F1 != +-F2, s = +-1, with
    /// [r]F1 + s [r]F2 = O
    pub cofactor_balance_c: Option<u64>,
}

/// x(P1 + P2) and x(P1 - P2) for P1 = (x1, y1), P2 = (x2, y2), x1 != x2, given 1/(x1 + x2).
#[inline]
fn pair_abscissae(f: &GF2n, a2: u64, x1: u64, y1: u64, x2: u64, y2: u64, inv: u64) -> (u64, u64) {
    let u = f.mul(y1 ^ y2, inv); // lambda of P1 + P2
    let v = f.mul(x2, inv); // lambda of P1 - P2 is u + v (y2 -> y2 + x2)
    let c = a2 ^ x1 ^ x2;
    let w = u ^ v;
    (f.sq(u) ^ u ^ c, f.sq(w) ^ w ^ c)
}

/// The exact 2-sum yield of the factor base: all unordered pairs of distinct abscissae, one
/// field inversion per batch of 4096 pairs (Montgomery's trick), abscissae collected, sorted and
/// deduplicated. With r (a prime dividing the order), the subgroup-restricted count and the
/// cofactor-balance count from h(P) = [r]P.
pub fn two_sum_yield(f: &GF2n, e: &BinCurve, fb: &FactorBase, r: Option<u128>) -> TwoSumYield {
    assert!(e.a4 == 0);
    let n = fb.xs.len();
    let pairs = (n as u64) * (n as u64 - 1) / 2;
    let mut xs: Vec<u64> = Vec::with_capacity((2 * pairs) as usize);
    const B: usize = 4096;
    let mut block: Vec<(usize, usize)> = Vec::with_capacity(B);
    let mut dens: Vec<u64> = Vec::with_capacity(B);
    let flush = |block: &mut Vec<(usize, usize)>, dens: &mut Vec<u64>, xs: &mut Vec<u64>| {
        let invs = f.batch_inv(dens);
        for (&(i, k), &inv) in block.iter().zip(&invs) {
            let (xp, xm) = pair_abscissae(f, e.a2, fb.xs[i], fb.ys[i], fb.xs[k], fb.ys[k], inv);
            xs.push(xp);
            xs.push(xm);
        }
        block.clear();
        dens.clear();
    };
    for i in 0..n {
        for k in i + 1..n {
            block.push((i, k));
            dens.push(fb.xs[i] ^ fb.xs[k]);
            if block.len() == B {
                flush(&mut block, &mut dens, &mut xs);
            }
        }
    }
    if !block.is_empty() {
        flush(&mut block, &mut dens, &mut xs);
    }
    xs.sort_unstable();
    xs.dedup();
    let count = xs.len() as u64;
    drop(xs);
    let mut out = TwoSumYield {
        pairs,
        count,
        ..Default::default()
    };
    let Some(r) = r else { return out };
    // h(P) = [r]P per abscissa; pairs with h(P1) = -+h(P2) share x(h)
    let hs: Vec<BPt> = (0..n).map(|i| e.mul(f, &fb.point(i), r)).collect();
    let mut groups: HashMap<(bool, u64), Vec<usize>> = HashMap::new();
    for (i, h) in hs.iter().enumerate() {
        let key = match *h {
            Pt::Inf => (true, 0),
            Pt::Aff(x, _) => (false, x),
        };
        groups.entry(key).or_default().push(i);
    }
    let mut sub_xs: Vec<u64> = vec![];
    let (mut combos, mut c) = (0u64, 0u64);
    for idx in groups.values() {
        for a in 0..idx.len() {
            for b in a + 1..idx.len() {
                let (i, k) = (idx[a], idx[b]);
                let plus = hs[i] == e.neg(&hs[k]); // P1 + P2 in the subgroup
                let minus = hs[i] == hs[k]; // P1 - P2 in the subgroup
                if !plus && !minus {
                    continue;
                }
                let inv = f.inv(fb.xs[i] ^ fb.xs[k]);
                let (xp, xm) = pair_abscissae(f, e.a2, fb.xs[i], fb.ys[i], fb.xs[k], fb.ys[k], inv);
                let m = u64::from(plus) + u64::from(minus);
                if plus {
                    sub_xs.push(xp);
                }
                if minus {
                    sub_xs.push(xm);
                }
                combos += m;
                if fb.xs[i] != 0 && fb.xs[k] != 0 {
                    // signs (e1, e2, s) with e1 h1 + s e2 h2 = 0: four per satisfied relation,
                    // times the two orders of (F1, F2)
                    c += 8 * m;
                }
            }
        }
    }
    sub_xs.sort_unstable();
    sub_xs.dedup();
    out.subgroup_count = Some(sub_xs.len() as u64);
    out.subgroup_combos = Some(combos);
    out.cofactor_balance_c = Some(c);
    out
}

/// Largest prime factor of the order and the cofactor.
pub fn largest_prime_factor(order: u128) -> Option<(u128, u128)> {
    if order > u64::MAX as u128 {
        return None;
    }
    let fac = factor_u64(order as u64);
    let r = fac.last()?.0 as u128;
    Some((r, order / r))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet, VecDeque};

    #[test]
    fn binary_volcano_frobenius_disc_decomposition() {
        // q = 2^12, t = 47: D = -14175 = -3^4 5^2 7, D_K = -7, f = 45
        let d = frobenius_disc(4096, 4097 - 47);
        assert_eq!((d.trace, d.d, d.d_k, d.f_pi), (47, -14175, -7, 45));
        assert_eq!(d.f_pi_factors, vec![(3, 2), (5, 1)]);
        assert_eq!(d.unfactored, None);
        assert_eq!(d.height(3), 2);
        assert_eq!(d.odd_heights(13), vec![(3, 2), (5, 1)]);
        // t = 115: D = -3159 = -3^5 13, D_K = -39, f = 9
        let d = frobenius_disc(4096, 4097 - 115);
        assert_eq!((d.d, d.d_k, d.f_pi), (-3159, -39, 9));
        // even conductor handling (generic q): q = 9, t = 0: D = -36 = 6^2 (-1), D_K = -4, f = 3
        let d = frobenius_disc(9, 10);
        assert_eq!((d.d, d.d_k, d.f_pi), (-36, -4, 3));
        assert_eq!(d.f_pi_factors, vec![(3, 1)]);
        // q = 7, t = 2: D = -24, D_K = -24, f = 1
        let d = frobenius_disc(7, 6);
        assert_eq!((d.d, d.d_k, d.f_pi), (-24, -24, 1));
        for &(dd, l, s) in &[
            (-7i128, 3u64, -1i32),
            (-39, 3, 0),
            (-7, 7, 0),
            (-3, 7, 1),
            (-39, 5, 1),
        ] {
            assert_eq!(kronecker_odd(dd, l), s, "({dd}/{l})");
        }
        // D = f^2 D_K with D_K fundamental and f = f_pi, for many orders
        for order in 4097 - 127..=4097 + 127 {
            let d = frobenius_disc(4096, order);
            assert_eq!(d.d, d.d_k * (d.f_pi * d.f_pi) as i128);
            assert!(d.d_k.rem_euclid(4) <= 1);
        }
        // only ell with (4q/ell) = 1 can divide D: at n = 51 that is ell = 7 alone
        assert_eq!(admissible_ells(1u128 << 51, 13), vec![7]);
        assert_eq!(admissible_ells(1u128 << 20, 13), vec![3, 5, 7, 11, 13]);
    }

    #[test]
    fn binary_volcano_roots_with_multiplicity() {
        let f = GF2n::new(13);
        let mut rng = Rng::new(1);
        let (a, b, c) = (5u64, 77u64, 1234u64);
        let p = poly::from_roots(&f, &[a, a, b, c, a]);
        let p = poly::scale(&f, &p, 99);
        assert_eq!(
            roots_with_multiplicity(&f, &p, &mut rng),
            vec![a, a, a, b, c]
        );
        // an irreducible quadratic factor contributes nothing
        let irr = vec![1u64, 1, 1]; // Y^2 + Y + 1 has no root when Tr(1) = 1 (n odd)
        let p = poly::mul(&f, &poly::from_roots(&f, &[b, b]), &irr);
        assert_eq!(roots_with_multiplicity(&f, &p, &mut rng), vec![b, b]);
    }

    #[test]
    fn binary_volcano_bsgs_order_matches_naive_at_n20() {
        let f = GF2n::new(20);
        let tau = trace_one_element(&f);
        assert_eq!(f.trace(tau), 1);
        let mut rng = Rng::new(20);
        let mut checked = 0;
        for _ in 0..4 {
            let a6 = loop {
                let v = f.random(&mut rng);
                if v != 0 {
                    break v;
                }
            };
            let mut orders = vec![];
            for a2 in [0u64, tau] {
                let e = BinCurve::new(a2, 0, a6);
                let naive = e.order_naive(&f);
                for seed in [1u64, 2] {
                    let mut r2 = Rng::new(seed);
                    assert_eq!(e.order(&f, &mut r2), naive, "a2 = {a2:#x}, a6 = {a6:#x}");
                    checked += 1;
                }
                orders.push(naive);
            }
            assert_eq!(orders[0] + orders[1], 2 * f.size() + 2, "twists");
            assert_ne!(orders[0], orders[1]);
            let d = frobenius_disc(f.size(), orders[0]);
            assert!(d.trace % 2 != 0, "ordinary: odd trace");
        }
        assert_eq!(checked, 16);
    }

    #[test]
    fn binary_volcano_same_order_a2_and_factor_base() {
        let f = GF2n::new(12);
        let tau = trace_one_element(&f);
        let mut rng = Rng::new(12);
        for _ in 0..5 {
            let j = f.random(&mut rng) | 1;
            let o0 = BinCurve::new(0, 0, f.inv(j)).order_naive(&f);
            let o1 = BinCurve::new(tau, 0, f.inv(j)).order_naive(&f);
            assert_eq!(same_order_a2(&f, j, o0, (0, tau), &mut rng), Some(0));
            assert_eq!(same_order_a2(&f, j, o1, (0, tau), &mut rng), Some(tau));
            // a wrong order is detected
            assert_eq!(same_order_a2(&f, j, o0 + 2, (0, tau), &mut rng), None);
            let e = BinCurve::new(tau, 0, f.inv(j));
            let fb = factor_base(&f, &e, 7);
            let brute: Vec<u64> = (0..128u64).filter(|&x| e.has_x(&f, x)).collect();
            assert_eq!(fb.xs, brute);
            assert_eq!(fb.points(), 2 * brute.len() - 1);
            for i in 0..fb.xs.len() {
                assert!(e.on_curve(&f, &fb.point(i)));
            }
        }
    }

    /// Brute force over all ordered pairs of factor-base points with the plain group law.
    fn naive_yield(f: &GF2n, e: &BinCurve, fb: &FactorBase, r: u128) -> TwoSumYield {
        let mut pts: Vec<BPt> = vec![];
        for i in 0..fb.xs.len() {
            let p = fb.point(i);
            pts.push(p);
            if fb.xs[i] != 0 {
                pts.push(e.neg(&p));
            }
        }
        let mut all = HashSet::new();
        let mut sub = HashSet::new();
        let (mut combos, mut c) = (0u64, 0u64);
        let hx = |p: &BPt| match *p {
            Pt::Inf => None,
            Pt::Aff(x, _) => Some(x),
        };
        for p1 in &pts {
            for p2 in &pts {
                if hx(p1) == hx(p2) {
                    continue; // P1 = +-P2
                }
                // both signs s: F1 + s F2 (the triples (F1, F2, s) of the definition)
                for p2s in [*p2, e.neg(p2)] {
                    let s = e.add(f, p1, &p2s);
                    all.insert(hx(&s).unwrap());
                    if e.mul(f, &s, r) == Pt::Inf {
                        sub.insert(hx(&s).unwrap());
                        let t0 = hx(p1) == Some(0) || hx(p2) == Some(0);
                        if !t0 {
                            c += 1;
                        }
                    }
                }
            }
        }
        // unordered abscissa pairs with the chosen representatives, both signs
        for i in 0..fb.xs.len() {
            for k in i + 1..fb.xs.len() {
                let (p1, p2) = (fb.point(i), fb.point(k));
                for q in [e.add(f, &p1, &p2), e.add(f, &p1, &e.neg(&p2))] {
                    if e.mul(f, &q, r) == Pt::Inf {
                        combos += 1;
                    }
                }
            }
        }
        let n = fb.xs.len() as u64;
        TwoSumYield {
            pairs: n * (n - 1) / 2,
            count: all.len() as u64,
            subgroup_count: Some(sub.len() as u64),
            subgroup_combos: Some(combos),
            cofactor_balance_c: Some(c),
        }
    }

    #[test]
    fn binary_volcano_two_sum_yield_matches_brute_force() {
        let mut rng = Rng::new(77);
        let mut nontrivial = 0;
        for n in [11u32, 12, 13] {
            let f = GF2n::new(n);
            let tau = trace_one_element(&f);
            for trial in 0..4 {
                let e = BinCurve::new(
                    if trial % 2 == 0 { 0 } else { tau },
                    0,
                    f.random(&mut rng) | 1,
                );
                let order = e.order_naive(&f);
                let (r, _) = largest_prime_factor(order).unwrap();
                let fb = factor_base(&f, &e, 5);
                let fast = two_sum_yield(&f, &e, &fb, Some(r));
                let slow = naive_yield(&f, &e, &fb, r);
                assert_eq!(fast, slow, "n = {n}, a6 = {:#x}", e.a6);
                assert!(fast.count <= 2 * fast.pairs);
                if fast.subgroup_count.unwrap() > 0 {
                    nontrivial += 1;
                }
                // without r the subgroup columns are absent and the rest unchanged
                let plain = two_sum_yield(&f, &e, &fb, None);
                assert_eq!((plain.pairs, plain.count), (fast.pairs, fast.count));
                assert_eq!(plain.subgroup_count, None);
            }
        }
        assert!(
            nontrivial > 0,
            "no curve with subgroup relations: test too weak"
        );
    }

    /// Every ordinary j over GF(2^12) with its two twists' orders; the full ell-graph of each
    /// class with v_ell(f_pi) >= 2 (ell odd <= 13) is built from Phi_ell mod 2 and the levels
    /// read off by breadth-first search from the floor are compared with the walk-based levels,
    /// the rational-isogeny counts with Kohel's description and with the kernel-polynomial count.
    #[test]
    fn binary_volcano_levels_match_brute_force_at_n12() {
        let f = GF2n::new(12);
        let q = f.size();
        let tau = trace_one_element(&f);
        let mut rng = Rng::new(12);
        // orders of (a2 = 0) for every j
        let mut order0: HashMap<u64, u128> = HashMap::new();
        let mut by_class: HashMap<u128, Vec<u64>> = HashMap::new(); // key: min(N, N')
        for j in 1..(1u64 << 12) {
            let o = BinCurve::new(0, 0, f.inv(j)).order_naive(&f);
            order0.insert(j, o);
            let o2 = 2 * q + 2 - o;
            by_class.entry(o.min(o2)).or_default().push(j);
        }
        let mut classes_checked = 0;
        let mut vertices_checked = 0;
        let mut craters = HashSet::new();
        for (&n_min, js) in &by_class {
            let disc = frobenius_disc(q, n_min);
            for (ell, h) in disc.odd_heights(13) {
                if h < 2 {
                    continue;
                }
                let g = EllGraph::new(&f, ell, h);
                let set: HashSet<u64> = js.iter().copied().collect();
                // full graph with multiplicity
                let mut nbrs: HashMap<u64, Vec<u64>> = HashMap::new();
                for &j in js {
                    let ns = g.neighbors(j, &mut rng);
                    for &x in &ns {
                        assert!(set.contains(&x), "neighbour leaves the class");
                    }
                    nbrs.insert(j, ns);
                }
                // brute-force levels: distance to the floor
                let mut dist: HashMap<u64, u32> = HashMap::new();
                let mut queue = VecDeque::new();
                for &j in js {
                    if nbrs[&j].len() == 1 {
                        dist.insert(j, 0);
                        queue.push_back(j);
                    }
                }
                while let Some(j) = queue.pop_front() {
                    let d = dist[&j];
                    for &x in &nbrs[&j] {
                        if !dist.contains_key(&x) {
                            dist.insert(x, d + 1);
                            queue.push_back(x);
                        }
                    }
                }
                let chi = kronecker_odd(disc.d_k, ell);
                for &j in js {
                    let d = dist[&j];
                    assert!(d <= h);
                    let level_bf = h - d;
                    assert_eq!(
                        g.level(j, &mut rng),
                        Some(level_bf),
                        "ell = {ell}, j = {j:#x}"
                    );
                    let cnt = nbrs[&j].len() as u64;
                    if level_bf == h {
                        assert_eq!(cnt, 1);
                    } else {
                        assert_eq!(cnt, ell + 1);
                    }
                    if level_bf == 0 {
                        let horizontal = nbrs[&j].iter().filter(|&&x| dist[&x] == h).count() as i32;
                        assert_eq!(horizontal, 1 + chi, "crater edges = 1 + (D_K/ell)");
                        craters.insert((n_min, ell, j));
                    }
                    vertices_checked += 1;
                }
                // kernel polynomials agree with the Phi multiplicity count on a few vertices
                for &j in js.iter().take(6) {
                    let e = BinCurve::new(0, 0, f.inv(j));
                    let ks = crate::binary::kernel_polys(&f, &e, ell, &mut rng);
                    assert_eq!(
                        ks.len(),
                        nbrs[&j].len(),
                        "kernel count, ell = {ell}, j = {j:#x}"
                    );
                    let mut targets: Vec<u64> = ks
                        .iter()
                        .map(|k| crate::binary::kohel_codomain(&f, &e, k).j(&f))
                        .collect();
                    targets.sort_unstable();
                    assert_eq!(targets, nbrs[&j]);
                }
                // the bounded component enumeration finds the whole volcano when not bounded
                let j0 = js[0];
                let comp = enumerate_component(
                    &g,
                    j0,
                    0,
                    order0[&j0],
                    (0, tau),
                    100_000,
                    &mut rng,
                    |_, _| {},
                );
                assert!(!comp.truncated);
                let mut seen: HashSet<u64> = HashSet::new();
                for v in &comp.vertices {
                    assert!(seen.insert(v.j));
                    assert!(dist.contains_key(&v.j));
                    assert_eq!(v.level, Some(h - dist[&v.j]));
                    assert_eq!(v.rational, nbrs[&v.j].len());
                    let a2 = v.a2.expect("twist found");
                    assert_eq!(
                        BinCurve::new(a2, 0, f.inv(v.j)).order_naive(&f),
                        order0[&j0]
                    );
                }
                // a connected component: everything reachable from j0 was found
                let mut reach = HashSet::from([j0]);
                let mut stack = vec![j0];
                while let Some(j) = stack.pop() {
                    for &x in &nbrs[&j] {
                        if reach.insert(x) {
                            stack.push(x);
                        }
                    }
                }
                assert_eq!(reach, seen);
                // and the bound is honoured
                let small =
                    enumerate_component(&g, j0, 0, order0[&j0], (0, tau), 3, &mut rng, |_, _| {});
                assert_eq!(small.vertices.len(), 3.min(reach.len()));
                assert_eq!(small.truncated, reach.len() > 3);
                classes_checked += 1;
            }
        }
        assert!(
            classes_checked >= 2,
            "classes with height >= 2: {classes_checked}"
        );
        assert!(vertices_checked >= 30, "vertices: {vertices_checked}");
        assert!(!craters.is_empty());
    }
}
