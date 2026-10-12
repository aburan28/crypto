//! Arm 3's branching rule: the minimum vertex cover of the monomial graph,
//! re-evaluated at every decision.
//!
//! **The graph.**  Vertices are the `3ℓ` coordinate bits of `X₁, X₂, X₃`; the
//! system is the library's descended symmetrised `S₄` (`weil_descend_s4`: the
//! correspondence equations `e_{i,d} = σ_{i,d}(x)` and the `n` Semaev
//! equations over the `e` variables).  Two vertices are adjacent iff some
//! monomial of degree ≥ 2 contains both.  A cover is a set of variables whose
//! assignment leaves every remaining monomial with at most one unassigned
//! variable, i.e. the residual system linear (the rule of WDSat's XG-ext).
//!
//! **Per node.**  The assigned bits are substituted (a monomial with an
//! assigned-0 bit vanishes, assigned-1 bits drop out, equal residual monomials
//! of one equation cancel), the graph of the unassigned bits is rebuilt from
//! the surviving monomials of degree ≥ 2, an exact minimum vertex cover is
//! computed (as the complement of a maximum independent set), and the next
//! decision is the cover vertex of largest residual degree (ties: lowest
//! lexicographic position), assigned 0 first.  Among the maximum independent
//! sets the one that is last in lexicographic order is used, so that the cover
//! prefers early variables.  When the residual graph has no edge the lowest
//! unassigned bit is taken.
//!
//! The `e` variables never share a monomial with an `x` bit in this system
//! (correspondence monomials are in `x`, Semaev monomials in `e`), so the
//! `e`-component of the graph is a separate component that no assignment of
//! `x` bits changes; it cannot influence which `x` bit is chosen and is only
//! measured (`graph-stats`), not covered per node.

use crypto_lib::cryptanalysis::binary_semaev_s4::{AnfF2m, AnfPoly, S4System};

/// Position of bit `j` of summand `s` in lexicographic order (summand-major,
/// most significant bit first): the SAT variable number minus one.
pub fn lex_pos(l: u32, s: usize, j: u32) -> usize {
    s * l as usize + (l - 1 - j) as usize
}

pub struct XGraph {
    pub nx: usize,
    pub l: u32,
    /// Per equation, the degree ≥ 2 monomials as masks over lexicographic positions.
    eqs: Vec<Vec<u64>>,
    pub n_monomials: usize,
    pub n_equations: usize,
}

/// Counts that characterise the graph of a whole system (target-independent for `x`).
#[derive(Debug, Clone)]
pub struct GraphShape {
    pub x_vertices: usize,
    pub x_edges: usize,
    pub x_edges_within_summand: usize,
    pub x_edges_across_summands: usize,
    pub x_min_vertex_cover: u32,
    pub e_vertices: usize,
    pub e_edges: usize,
}

impl XGraph {
    /// The correspondence equations `σ_{i,d}(x)` exactly as `weil_descend_s4` builds them (its half 1),
    /// without the Semaev half, which the branching rule never reads.
    pub fn new(l: u32) -> Self {
        let xs: Vec<AnfF2m> = (0..3)
            .map(|i| AnfF2m::from_vars(i as u32 * l, l as usize))
            .collect();
        let sigma1 = xs[0].xor(&xs[1]).xor(&xs[2]);
        let x0x1 = xs[0].mul(&xs[1]);
        let x0x2 = xs[0].mul(&xs[2]);
        let x1x2 = xs[1].mul(&xs[2]);
        let sigma2 = x0x1.xor(&x0x2).xor(&x1x2);
        let sigma3 = x0x1.mul(&xs[2]);
        let rows: Vec<Vec<AnfPoly>> = [sigma1, sigma2, sigma3]
            .iter()
            .map(|s| s.coeffs.clone())
            .collect();
        Self::from_rows(l, &rows)
    }

    pub fn from_system(sys: &S4System) -> Self {
        Self::from_rows(sys.l, &sys.correspondence)
    }

    fn from_rows(l: u32, rows: &[Vec<AnfPoly>]) -> Self {
        let nx = 3 * l as usize;
        assert!(nx <= 63, "monomial masks are u64");
        let mut eqs = Vec::new();
        let mut n_monomials = 0usize;
        for fam in rows {
            for poly in fam {
                let mut masks = Vec::new();
                for mono in poly.monomials() {
                    if mono.len() >= 2 {
                        let mut m = 0u64;
                        for &v in mono {
                            let (s, j) = (v / l, v % l);
                            m |= 1u64 << lex_pos(l, s as usize, j);
                        }
                        masks.push(m);
                    }
                }
                n_monomials += masks.len();
                if !masks.is_empty() {
                    eqs.push(masks);
                }
            }
        }
        let n_equations = eqs.len();
        Self {
            nx,
            l,
            eqs,
            n_monomials,
            n_equations,
        }
    }

    /// Adjacency of the unassigned bits after substituting `assigned` (mask) with values `ones`.
    pub fn adjacency(&self, assigned: u64, ones: u64, adj: &mut [u64], buf: &mut Vec<u64>) {
        for a in adj.iter_mut().take(self.nx) {
            *a = 0;
        }
        let zeros = assigned & !ones;
        for eq in &self.eqs {
            buf.clear();
            for &m in eq {
                if m & zeros != 0 {
                    continue;
                }
                let r = m & !assigned;
                if r.count_ones() >= 2 {
                    buf.push(r);
                }
            }
            buf.sort_unstable();
            let mut i = 0;
            while i < buf.len() {
                let mut j = i;
                while j < buf.len() && buf[j] == buf[i] {
                    j += 1;
                }
                if (j - i) % 2 == 1 {
                    let r = buf[i];
                    let mut bits = r;
                    while bits != 0 {
                        let v = bits.trailing_zeros() as usize;
                        bits &= bits - 1;
                        adj[v] |= r & !(1u64 << v);
                    }
                }
                i = j;
            }
        }
    }
}

/// Exact size of a maximum independent set of the subgraph on `avail`.
pub fn alpha(avail: u64, adj: &[u64], steps: &mut u64) -> u32 {
    *steps += 1;
    if avail == 0 {
        return 0;
    }
    let mut best_v = usize::MAX;
    let mut best_d = 0u32;
    let mut bits = avail;
    while bits != 0 {
        let v = bits.trailing_zeros() as usize;
        bits &= bits - 1;
        let d = (adj[v] & avail).count_ones();
        if d > best_d {
            best_d = d;
            best_v = v;
        }
    }
    if best_d == 0 {
        return avail.count_ones();
    }
    let v = best_v;
    let without = alpha(avail & !(1u64 << v), adj, steps);
    let with = 1 + alpha(avail & !(adj[v] | (1u64 << v)), adj, steps);
    without.max(with)
}

/// The maximum independent set that is last in lexicographic order (greedy
/// from the highest position down, each vertex kept iff an optimum still
/// contains it together with the vertices kept so far).
pub fn lex_last_mis(vertices: u64, adj: &[u64], steps: &mut u64) -> u64 {
    let target = alpha(vertices, adj, steps);
    let (mut chosen, mut avail) = (0u64, vertices);
    let mut bits = vertices;
    while bits != 0 {
        let v = 63 - bits.leading_zeros() as usize;
        bits &= !(1u64 << v);
        if avail & (1u64 << v) == 0 {
            continue;
        }
        let rest = avail & !(adj[v] | (1u64 << v));
        if chosen.count_ones() + 1 + alpha(rest, adj, steps) == target {
            chosen |= 1u64 << v;
            avail = rest;
        } else {
            avail &= !(1u64 << v);
        }
    }
    chosen
}

#[derive(Default, Clone, Debug)]
pub struct DynStats {
    pub decisions: u64,
    /// Decisions whose variable differs from the lexicographically first unassigned bit.
    pub differs_from_lex: u64,
    pub no_edge_fallbacks: u64,
    pub cover_size_sum: u64,
    pub mis_steps: u64,
    pub max_mis_steps_one_node: u64,
}

pub struct MvcDecider<'g> {
    g: &'g XGraph,
    adj: Vec<u64>,
    buf: Vec<u64>,
    pub stats: DynStats,
}

impl<'g> MvcDecider<'g> {
    pub fn new(g: &'g XGraph) -> Self {
        Self {
            g,
            adj: vec![0; g.nx],
            buf: Vec::new(),
            stats: DynStats::default(),
        }
    }

    /// Lexicographic position of the next bit to decide, or `None` when all bits are assigned.
    pub fn choose(&mut self, asg: &[Option<bool>]) -> Option<usize> {
        let nx = self.g.nx;
        let (mut assigned, mut ones) = (0u64, 0u64);
        for (k, a) in asg.iter().take(nx).enumerate() {
            if let Some(v) = a {
                assigned |= 1u64 << k;
                if *v {
                    ones |= 1u64 << k;
                }
            }
        }
        let full = (1u64 << nx) - 1;
        let unassigned = full & !assigned;
        if unassigned == 0 {
            return None;
        }
        let lex_first = unassigned.trailing_zeros() as usize;
        self.stats.decisions += 1;
        self.g
            .adjacency(assigned, ones, &mut self.adj, &mut self.buf);
        let mut active = 0u64;
        let mut bits = unassigned;
        while bits != 0 {
            let v = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            if self.adj[v] & unassigned != 0 {
                active |= 1u64 << v;
            }
        }
        let pick = if active == 0 {
            self.stats.no_edge_fallbacks += 1;
            lex_first
        } else {
            let mut steps = 0u64;
            let mis = lex_last_mis(active, &self.adj, &mut steps);
            self.stats.mis_steps += steps;
            self.stats.max_mis_steps_one_node = self.stats.max_mis_steps_one_node.max(steps);
            let cover = active & !mis;
            self.stats.cover_size_sum += cover.count_ones() as u64;
            let (mut best, mut best_d) = (usize::MAX, 0u32);
            let mut bits = cover;
            while bits != 0 {
                let v = bits.trailing_zeros() as usize;
                bits &= bits - 1;
                let d = (self.adj[v] & unassigned).count_ones();
                if best == usize::MAX || d > best_d {
                    best = v;
                    best_d = d;
                }
            }
            best
        };
        if pick != lex_first {
            self.stats.differs_from_lex += 1;
        }
        Some(pick)
    }
}

/// Shape of the monomial graph of a whole system (diagnostic for `graph-stats`).
pub fn shape(sys: &S4System) -> GraphShape {
    let g = XGraph::from_system(sys);
    let l = sys.l as usize;
    let mut adj = vec![0u64; g.nx];
    let mut buf = Vec::new();
    g.adjacency(0, 0, &mut adj, &mut buf);
    let (mut within, mut across) = (0usize, 0usize);
    for v in 0..g.nx {
        let mut bits = (adj[v] >> (v + 1)) << (v + 1);
        while bits != 0 {
            let w = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            if v / l == w / l {
                within += 1;
            } else {
                across += 1;
            }
        }
    }
    let all = (1u64 << g.nx) - 1;
    let mut steps = 0u64;
    let alpha_x = alpha(all, &adj, &mut steps);
    // The e-component: the Semaev equations are monomials in e only.
    let n_e = sys.n_e_vars() as usize;
    assert!(n_e <= 64, "e-graph masks are u64");
    let mut eadj = vec![0u64; n_e];
    for p in &sys.semaev {
        for mono in p.monomials() {
            if mono.len() >= 2 {
                for &a in mono {
                    for &b in mono {
                        if a != b {
                            eadj[a as usize] |= 1u64 << b;
                        }
                    }
                }
            }
        }
    }
    let e_edges = eadj.iter().map(|m| m.count_ones() as usize).sum::<usize>() / 2;
    GraphShape {
        x_vertices: g.nx,
        x_edges: within + across,
        x_edges_within_summand: within,
        x_edges_across_summands: across,
        x_min_vertex_cover: g.nx as u32 - alpha_x,
        e_vertices: n_e,
        e_edges,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn complete_multipartite(parts: &[usize]) -> (u64, Vec<u64>) {
        let n: usize = parts.iter().sum();
        let mut adj = vec![0u64; n];
        let mut part_of = Vec::new();
        for (p, &c) in parts.iter().enumerate() {
            for _ in 0..c {
                part_of.push(p);
            }
        }
        for v in 0..n {
            for w in 0..n {
                if part_of[v] != part_of[w] {
                    adj[v] |= 1u64 << w;
                }
            }
        }
        ((1u64 << n) - 1, adj)
    }

    #[test]
    fn alpha_of_complete_multipartite_is_the_largest_part() {
        let (all, adj) = complete_multipartite(&[3, 5, 4]);
        let mut s = 0;
        assert_eq!(alpha(all, &adj, &mut s), 5);
        // Last in lexicographic order: the part at the highest positions among the optima.
        let mis = lex_last_mis(all, &adj, &mut s);
        assert_eq!(mis.count_ones(), 5);
        assert_eq!(mis, 0b11111 << 3);
        let (all, adj) = complete_multipartite(&[4, 4, 4]);
        let mis = lex_last_mis(all, &adj, &mut s);
        assert_eq!(mis, 0b1111 << 8, "ties go to the last part");
    }

    #[test]
    fn alpha_agrees_with_brute_force_on_random_graphs() {
        let mut rng = 0x1234_5678_9abc_def0u64;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        for _ in 0..200 {
            let n = 4 + (next() % 9) as usize;
            let mut adj = vec![0u64; n];
            for v in 0..n {
                for w in v + 1..n {
                    if next() % 3 == 0 {
                        adj[v] |= 1 << w;
                        adj[w] |= 1 << v;
                    }
                }
            }
            let all = (1u64 << n) - 1;
            let mut brute = 0;
            for s in 0..=all {
                if (0..n).all(|v| s >> v & 1 == 0 || adj[v] & s == 0) {
                    brute = brute.max(s.count_ones());
                }
            }
            let mut st = 0;
            assert_eq!(alpha(all, &adj, &mut st), brute);
            let mis = lex_last_mis(all, &adj, &mut st);
            assert_eq!(mis.count_ones(), brute);
            assert!((0..n).all(|v| mis >> v & 1 == 0 || adj[v] & mis == 0));
        }
    }
}
