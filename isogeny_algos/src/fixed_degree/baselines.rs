//! Classical baselines for "find an isogeny E1 -> E2 of a given degree d" between supersingular
//! curves over F_{p^2} (u64 p, p = 3 mod 4), the problem of Benčina–Kutas–Merz–Petit–Stopar–
//! Weitkämper 2023, against which their algorithm is compared.
//!
//! * [`exhaustive_prime_degree_phi`], [`exhaustive_prime_degree_kernel`]: prime d, enumerate all
//!   d + 1 isogenies of degree d out of E1 (roots of Phi_d(j(E1), Y), or the d + 1 cyclic
//!   subgroups of E1[d] pushed through Vélu) and compare with j(E2).
//! * [`mitm_smooth`]: d = l^e, meet in the middle on j-paths: all non-backtracking l-paths of
//!   length ceil(e/2) from j1 and floor(e/2) from j2, matched by a hash map.
//! * [`vow_claw`]: van Oorschot–Wiener parallel collision search (run serially) for the claw
//!   f1(x) = f2(y) between the two half-path functions, with distinguished points and function
//!   versions.
//! * [`grover_tani_cost_model`]: query counts of Grover and Tani's quantum claw finding; *not* a
//!   simulation.
//!
//! **Paths, cyclicity and j = 0, 1728.** All searches work on the j-line with the *distinct*
//! roots of Phi_l(j, Y) (sorted, which gives the canonical neighbour order used by vOW). A
//! j-path j_0 .. j_e with Phi_l(j_i, j_{i+1}) = 0 and j_{i+1} != j_{i-1} is the j-path of a cyclic
//! isogeny of degree l^e: each step is an l-isogeny phi_i, and phi_{i+1} cannot be (an isomorphism
//! composed with) the dual of phi_i because the dual lands on j_{i-1}. The converse is false at
//! multiple edges: if Phi_l(j, Y) has a repeated root (two non-equivalent l-isogenies to the same
//! j', e.g. out of j = 0 or 1728 where Aut(E) is larger, or at genuine double edges), cyclic walks
//! that return to j_{i-1} along the *other* edge, and the multiplicity of edges, are not
//! represented. So the search space is the set of j-non-backtracking walks, a subset (modulo
//! automorphisms) of the cyclic l^e-isogenies; every path returned is cyclic. For p >> l the
//! affected vertices are a vanishing fraction (j = 0, 1728 and O(l) double-edge vertices).
//! At a vertex with fewer than the generic number of non-backtracking neighbours the vOW index
//! digit is reduced modulo the actual count (non-uniform there) and a vertex with none ends the
//! path as a "dead" output (it only produces useless collisions).
use crate::curve::*;
use crate::field::{Field, Rng, Zp2};
use crate::find::modpoly::Phi;
use crate::path::graph::Path;
use std::collections::HashMap;

type J = (u64, u64);

// ------------------------------------------------------------------ neighbour cache

/// Sorted distinct roots of Phi_l(j, Y) over F_{p^2}, memoised. `root_findings` counts the
/// distinct j for which roots were actually computed (the memo only saves wall time; the cost
/// models below count graph steps / function evaluations, not cache misses).
pub struct Neighbours<'a> {
    pub f2: &'a Zp2,
    pub phi: &'a Phi<Zp2>,
    pub map: HashMap<J, Vec<J>>,
    pub root_findings: usize,
    rng: Rng,
}
impl<'a> Neighbours<'a> {
    pub fn new(f2: &'a Zp2, phi: &'a Phi<Zp2>) -> Self {
        Neighbours {
            f2,
            phi,
            map: HashMap::new(),
            root_findings: 0,
            rng: Rng::new(0x5eed),
        }
    }
    pub fn get(&mut self, j: J) -> &Vec<J> {
        if !self.map.contains_key(&j) {
            let mut v = self.phi.neighbors(self.f2, j, &mut self.rng);
            v.sort();
            v.dedup();
            self.root_findings += 1;
            self.map.insert(j, v);
        }
        &self.map[&j]
    }
    /// Non-backtracking neighbours of j coming from `prev`, in canonical (sorted) order.
    pub fn forward(&mut self, j: J, prev: Option<J>) -> Vec<J> {
        self.get(j)
            .iter()
            .copied()
            .filter(|&n| Some(n) != prev)
            .collect()
    }
}

/// Phi_l on every edge and j_{i+1} != j_{i-1}: the j-path of a cyclic l^len isogeny (see the
/// module documentation).
pub fn verify_cyclic_path(f2: &Zp2, phi: &Phi<Zp2>, path: &Path<J>) -> bool {
    let n = path.js.len();
    if path.ells.len() + 1 != n || path.ells.iter().any(|&l| l != phi.ell) {
        return false;
    }
    for i in 0..n.saturating_sub(1) {
        if phi.eval(f2, path.js[i], path.js[i + 1]) != f2.zero() {
            return false;
        }
        if i >= 1 && path.js[i + 1] == path.js[i - 1] {
            return false;
        }
    }
    true
}

// ------------------------------------------------------------------ exhaustive, prime degree

/// Result of [`exhaustive_prime_degree_phi`].
#[derive(Clone, Debug)]
pub struct PhiExhaustive {
    /// Distinct roots of Phi_d(j1, Y) (d + 1 counted with multiplicity).
    pub candidates: usize,
    /// Position of j2 among the sorted roots.
    pub match_index: Option<usize>,
}

/// Is there a d-isogeny j1 -> j2 (d prime)? Roots of Phi_d(j1, Y) over F_{p^2}, compared with j2.
/// Cost: one root finding of a degree-(d+1) polynomial, i.e. O~(d) field operations after Phi_d.
pub fn exhaustive_prime_degree_phi(f2: &Zp2, phi: &Phi<Zp2>, j1: J, j2: J) -> PhiExhaustive {
    let mut nb = Neighbours::new(f2, phi);
    let roots = nb.get(j1).clone();
    PhiExhaustive {
        candidates: roots.len(),
        match_index: roots.iter().position(|&r| r == j2),
    }
}

/// A d-isogeny E1 -> E2 found by [`exhaustive_prime_degree_kernel`].
#[derive(Clone, Debug)]
pub struct KernelMatch {
    /// Generator of the kernel (a point of order d on E1 over F_{p^2}).
    pub kernel_gen: Pt<J>,
    /// Vélu codomain (isomorphic to E2: same j).
    pub codomain: Curve<J>,
    /// Kernels tried before (and including) the match; at most d + 1.
    pub tried: usize,
}

/// Enumerate the d + 1 cyclic subgroups of order d of E1[d] (d an odd prime dividing p + 1, E1
/// supersingular over F_{p^2} with E1(F_{p^2}) = (Z/(p+1))^2, e.g. any curve reached from
/// y^2 = x^3 + x by F_{p^2}-rational Vélu isogenies) and return the first whose Vélu codomain has
/// j = j2. Returns None if no subgroup matches, or if E1(F_{p^2}) is not killed by p + 1 (the
/// quadratic twist, group (Z/(p-1))^2: use a twist or d | p - 1 instead; not handled).
/// Cost: (d + 1) Vélu isogenies of degree d, O(d^2) field operations.
pub fn exhaustive_prime_degree_kernel(
    f2: &Zp2,
    e1: &Curve<J>,
    p: u64,
    d: u64,
    j2: J,
    rng: &mut Rng,
) -> Option<KernelMatch> {
    use crate::kernel::velu::velu_from_point;
    assert!(
        d % 2 == 1 && (p + 1).is_multiple_of(d),
        "need odd d | p + 1"
    );
    let r = random_point_f(f2, e1, rng);
    if pmul(f2, e1, &r, (p + 1) as u128) != Pt::Inf {
        return None;
    }
    let pt = crate::testdata::ss_point_of_order(f2, e1, p, d, rng);
    let mut mult = vec![Pt::Inf];
    for _ in 1..d {
        let last = *mult.last().unwrap();
        mult.push(padd(f2, e1, &last, &pt));
    }
    let q = loop {
        let q = crate::testdata::ss_point_of_order(f2, e1, p, d, rng);
        if !mult.contains(&q) {
            break q;
        }
    };
    // subgroups <Q>, <P + i Q> for i = 0..d-1
    let mut gens = vec![q];
    let mut cur = pt;
    for _ in 0..d {
        gens.push(cur);
        cur = padd(f2, e1, &cur, &q);
    }
    for (t, g) in gens.iter().enumerate() {
        let iso = velu_from_point(f2, e1, g, d);
        if jinv(f2, &iso.cod) == j2 {
            return Some(KernelMatch {
                kernel_gen: *g,
                codomain: iso.cod,
                tried: t + 1,
            });
        }
    }
    None
}

// ------------------------------------------------------------------ meet in the middle

/// Counters for [`mitm_smooth`].
#[derive(Clone, Debug, Default)]
pub struct MitmStats {
    /// Vertices expanded (one neighbour computation each) on the two sides.
    pub expansions: usize,
    /// Stored leaves (side 1) and probed leaves (side 2).
    pub leaves1: usize,
    pub leaves2: usize,
    /// Distinct j for which Phi_l roots were computed.
    pub root_findings: usize,
}

#[derive(Clone, Copy)]
struct Node {
    j: J,
    prev: Option<J>,
    parent: usize,
}

fn tree(nb: &mut Neighbours, root: J, depth: usize, st: &mut MitmStats) -> Vec<Vec<Node>> {
    let mut levels = vec![vec![Node {
        j: root,
        prev: None,
        parent: usize::MAX,
    }]];
    for _ in 0..depth {
        let mut next = vec![];
        for (i, n) in levels.last().unwrap().iter().enumerate() {
            st.expansions += 1;
            for c in nb.forward(n.j, n.prev) {
                next.push(Node {
                    j: c,
                    prev: Some(n.j),
                    parent: i,
                });
            }
        }
        levels.push(next);
    }
    levels
}

fn trace(levels: &[Vec<Node>], mut idx: usize) -> Vec<J> {
    let mut out = vec![];
    for lv in (0..levels.len()).rev() {
        let n = levels[lv][idx];
        out.push(n.j);
        idx = n.parent;
    }
    out.reverse();
    out
}

/// Meet in the middle for a cyclic l^e-isogeny j1 -> j2: all non-backtracking j-paths of length
/// ceil(e/2) from j1 (stored by endpoint) and floor(e/2) from j2 (probed); a match whose two
/// last steps differ (so the joined path does not backtrack at the meeting vertex) gives a
/// non-backtracking path of length e. Work and memory ~ (l+1) l^(ceil(e/2)-1) vertices per side.
pub fn mitm_smooth(
    f2: &Zp2,
    phi: &Phi<Zp2>,
    e: usize,
    j1: J,
    j2: J,
) -> (Option<Path<J>>, MitmStats) {
    let ell = phi.ell;
    let mut st = MitmStats::default();
    let mut nb = Neighbours::new(f2, phi);
    let (k1, k2) = (e.div_ceil(2), e / 2);
    let t1 = tree(&mut nb, j1, k1, &mut st);
    let leaves1 = t1.last().unwrap();
    st.leaves1 = leaves1.len();
    let mut table: HashMap<J, Vec<usize>> = HashMap::new();
    for (i, n) in leaves1.iter().enumerate() {
        table.entry(n.j).or_default().push(i);
    }
    let t2 = tree(&mut nb, j2, k2, &mut st);
    let leaves2 = t2.last().unwrap();
    st.leaves2 = leaves2.len();
    st.root_findings = nb.root_findings;
    for (i2, n2) in leaves2.iter().enumerate() {
        let Some(c) = table.get(&n2.j) else { continue };
        for &i1 in c {
            let n1 = leaves1[i1];
            if n1.prev.is_some() && n2.prev.is_some() && n1.prev == n2.prev {
                continue;
            }
            let mut js = trace(&t1, i1);
            let back = trace(&t2, i2);
            js.extend(back.iter().rev().skip(1));
            let path = Path {
                ells: vec![ell; js.len() - 1],
                js,
            };
            return (Some(path), st);
        }
    }
    (None, st)
}

// ------------------------------------------------------------------ van Oorschot-Wiener

/// Parameters for [`vow_claw`].
#[derive(Clone, Debug)]
pub struct VowParams {
    /// Memory w: number of distinguished-point slots.
    pub memory: usize,
    /// Distinguished points per function version (vOW recommend 10 w).
    pub dist_per_version: usize,
    /// Distinguished-point fraction theta = alpha sqrt(w / |S|) (vOW: alpha = 2.25), capped at 1.
    pub alpha: f64,
    /// Give up after this many function evaluations.
    pub max_evals: u64,
    pub seed: u64,
}
impl VowParams {
    pub fn new(memory: usize) -> Self {
        VowParams {
            memory,
            dist_per_version: 10 * memory,
            alpha: 2.25,
            max_evals: 1 << 26,
            seed: 1,
        }
    }
}

/// Counters for [`vow_claw`].
#[derive(Clone, Debug, Default)]
pub struct VowStats {
    /// Evaluations of f (each walks one half-path of up to ceil(e/2) steps).
    pub f_evals: u64,
    pub distinguished: u64,
    /// Collisions located (two different inputs with the same g output).
    pub collisions: u64,
    /// Collisions that were not the wanted claw (same side, hash collision, or backtracking at
    /// the meeting vertex).
    pub useless: u64,
    /// Trails abandoned for exceeding 20 / theta steps.
    pub aborted_trails: u64,
    /// Function versions used (1 + number of restarts).
    pub versions: u64,
    /// |S| = 2 N (two sides times N indices).
    pub space: u64,
    pub theta: f64,
    pub root_findings: usize,
}

fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

struct Claw<'a, 'b> {
    nb: &'b mut Neighbours<'a>,
    ell: u64,
    k: [usize; 2],
    root: [J; 2],
    n: [u64; 2],
    space: u64,
    evals: u64,
}

/// The decoded half-path: the j's and whether it ended early (no forward neighbour).
struct Half {
    js: Vec<J>,
    dead: bool,
}

impl Claw<'_, '_> {
    /// Path from root[side] encoded by idx: first digit base l+1, then base l; each digit picks
    /// (modulo the actual count) among the sorted non-backtracking neighbours.
    fn half(&mut self, side: usize, idx: u64) -> Half {
        let mut idx = idx % self.n[side];
        let mut js = vec![self.root[side]];
        let mut prev = None;
        for s in 0..self.k[side] {
            let base = if s == 0 { self.ell + 1 } else { self.ell };
            let digit = idx % base;
            idx /= base;
            let cur = *js.last().unwrap();
            let fw = self.nb.forward(cur, prev);
            if fw.is_empty() {
                return Half { js, dead: true };
            }
            prev = Some(cur);
            js.push(fw[(digit % fw.len() as u64) as usize]);
        }
        Half { js, dead: false }
    }
    fn f(&mut self, x: u64) -> Half {
        self.evals += 1;
        let (side, idx) = ((x >= self.n[0]) as usize, x % self.n[0]);
        self.half(side, idx)
    }
    fn g(&mut self, x: u64, v: u64) -> u64 {
        let h = self.f(x);
        let j = *h.js.last().unwrap();
        let tag = if h.dead { x ^ 0xdead_beef_0000 } else { 0 };
        mix(mix(j.0 ^ v.wrapping_mul(0x1234_5678_9abc_def1)) ^ j.1 ^ tag.rotate_left(17))
            % self.space
    }
}

/// van Oorschot–Wiener claw search for a cyclic l^e-isogeny j1 -> j2 (see the module notes).
/// S = {0, 1} x [0, N) with N = (l+1) l^(k1-1), k1 = ceil(e/2); f(0, x) = endpoint of the
/// half-path from j1 encoded by x, f(1, x) = endpoint of the half-path of length k2 = floor(e/2)
/// from j2 encoded by x mod N2 (N2 = (l+1) l^(k2-1); for odd e this makes each side-2 path appear
/// l times). g_v(x) = h_v(f(x)) with a hash h_v onto S that depends on the function version v.
/// Trails run until a distinguished point (hash below theta) and are stored in a w-slot table;
/// two trails meeting in one distinguished point are walked to their collision. A collision
/// f(0, x) = f(1, y) whose half-paths do not share the last vertex is the claw, i.e. a solution.
/// The version changes after `dist_per_version` distinguished points (and the table is cleared).
/// Expected work for a single claw: about 2.5 sqrt(|S|^3 / w) evaluations of f.
pub fn vow_claw(
    f2: &Zp2,
    phi: &Phi<Zp2>,
    e: usize,
    j1: J,
    j2: J,
    prm: &VowParams,
) -> (Option<Path<J>>, VowStats) {
    let ell = phi.ell as u64;
    let mut st = VowStats::default();
    let (k1, k2) = (e.div_ceil(2), e / 2);
    let size = |k: usize| {
        if k == 0 {
            1
        } else {
            (ell + 1) * ell.pow(k as u32 - 1)
        }
    };
    let n = [size(k1), size(k2)];
    let space = 2 * n[0];
    st.space = space;
    let mut nb = Neighbours::new(f2, phi);
    if e == 0 {
        return (
            (j1 == j2).then(|| Path {
                js: vec![j1],
                ells: vec![],
            }),
            st,
        );
    }
    let mut c = Claw {
        nb: &mut nb,
        ell,
        k: [k1, k2],
        root: [j1, j2],
        n,
        space,
        evals: 0,
    };
    let w = prm.memory.max(1);
    let theta = (prm.alpha * (w as f64 / space as f64).sqrt()).min(1.0);
    st.theta = theta;
    let thr = if theta >= 1.0 {
        u64::MAX
    } else {
        (theta * u64::MAX as f64) as u64
    };
    let max_trail = (20.0 / theta).ceil() as u64 + 1;
    let mut rng = Rng::new(prm.seed);
    let mut v = 0u64;
    st.versions = 1;
    let mut table: Vec<Option<(u64, u64, u64)>> = vec![None; w];
    let mut dist_this_version = 0usize;
    let is_dist = |x: u64, v: u64| mix(x ^ mix(v ^ 0x00d1_57d1)) <= thr;
    let mut result = None;
    'outer: while c.evals < prm.max_evals {
        if dist_this_version >= prm.dist_per_version {
            v += 1;
            st.versions += 1;
            dist_this_version = 0;
            table.iter_mut().for_each(|s| *s = None);
        }
        let start = rng.below(space);
        let mut x = start;
        let mut len = 0u64;
        let mut ok = false;
        while len < max_trail {
            if is_dist(x, v) && len > 0 {
                ok = true;
                break;
            }
            x = c.g(x, v);
            len += 1;
        }
        if !ok {
            st.aborted_trails += 1;
            continue;
        }
        st.distinguished += 1;
        dist_this_version += 1;
        let slot = (mix(x ^ v.rotate_left(32)) % w as u64) as usize;
        let prev = table[slot];
        table[slot] = Some((start, x, len));
        let Some((s2, d2, l2)) = prev else { continue };
        if d2 != x || s2 == start {
            continue;
        }
        // locate the collision: align the trails, then step together
        let (mut a, mut la, mut b, mut lb) = (start, len, s2, l2);
        while la > lb {
            a = c.g(a, v);
            la -= 1;
        }
        while lb > la {
            b = c.g(b, v);
            lb -= 1;
        }
        if a == b {
            continue; // one start lies on the other trail (Robin Hood)
        }
        loop {
            let (ga, gb) = (c.g(a, v), c.g(b, v));
            if ga == gb {
                break;
            }
            a = ga;
            b = gb;
        }
        st.collisions += 1;
        let (ha, hb) = (c.f(a), c.f(b));
        let (sa, sb) = ((a >= n[0]) as usize, (b >= n[0]) as usize);
        let ja = *ha.js.last().unwrap();
        let jb = *hb.js.last().unwrap();
        if sa == sb || ha.dead || hb.dead || ja != jb {
            st.useless += 1;
            continue;
        }
        let (h1, h2) = if sa == 0 { (ha, hb) } else { (hb, ha) };
        let p1 = (h1.js.len() >= 2).then(|| h1.js[h1.js.len() - 2]);
        let p2 = (h2.js.len() >= 2).then(|| h2.js[h2.js.len() - 2]);
        if p1.is_some() && p1 == p2 {
            st.useless += 1;
            continue;
        }
        let mut js = h1.js.clone();
        js.extend(h2.js.iter().rev().skip(1));
        result = Some(Path {
            ells: vec![ell as usize; js.len() - 1],
            js,
        });
        break 'outer;
    }
    st.f_evals = c.evals;
    st.root_findings = nb.root_findings;
    (result, st)
}

/// The vOW estimate 2.5 sqrt(|S|^3 / w) of f-evaluations to find one golden collision.
pub fn vow_expected_evals(space: u64, memory: usize) -> f64 {
    2.5 * ((space as f64).powi(3) / memory as f64).sqrt()
}

// ------------------------------------------------------------------ quantum cost model

/// Oracle-query counts for quantum search (cost model only; nothing is simulated).
#[derive(Clone, Debug)]
pub struct QuantumCost {
    /// Grover over all candidates: (pi / 4) sqrt(#candidates) queries (one marked element).
    pub grover_queries: f64,
    /// Tani's claw finding between functions on domains of sizes N1, N2:
    /// (N1 N2)^(1/3) queries (for N1 <= N2 <= N1^2; constant factors ignored).
    pub tani_queries: f64,
}

/// Query counts used in reports for d = l^e: the candidates are the (l+1) l^(e-1) cyclic
/// l^e-isogenies out of E1, the claw domains the half-paths of length ceil(e/2) and floor(e/2).
/// Only oracle queries are counted (no circuit, memory or QRAM cost, no error correction).
pub fn grover_tani_cost_model(ell: u64, e: u32) -> QuantumCost {
    let l = ell as f64;
    let size = |k: u32| {
        if k == 0 {
            1.0
        } else {
            (l + 1.0) * l.powi(k as i32 - 1)
        }
    };
    let cand = size(e);
    let (n1, n2) = (size(e.div_ceil(2)), size(e / 2));
    QuantumCost {
        grover_queries: std::f64::consts::FRAC_PI_4 * cand.sqrt(),
        tani_queries: (n1 * n2).cbrt(),
    }
}

/// The same model for an explicit candidate count and claw domain sizes.
pub fn grover_tani_cost(candidates: f64, n1: f64, n2: f64) -> QuantumCost {
    QuantumCost {
        grover_queries: std::f64::consts::FRAC_PI_4 * candidates.sqrt(),
        tani_queries: (n1 * n2).cbrt(),
    }
}
