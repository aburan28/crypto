//! `MultiSolve` with `OracleT`: every node calls `GroebnerSafe`; wild nodes
//! branch on the next summand coordinate, with the residual-weight bound.

use crate::boolpoly::{Poly, W};
use crate::f4::{linear_solutions, truncated_groebner, Outcome};

#[derive(Clone, Copy)]
pub struct Config {
    pub d_max: u32,
    pub budget_xor: u64,
    pub max_calls: u64,
    pub max_free: u32,
    pub matrix_cap_words: u64,
}

pub struct Instance {
    pub eqs: Vec<Poly>,
    pub n_vars: usize,
    /// Branching variables in order, with their summand index.
    pub branch: Vec<(usize, usize)>,
    /// Coordinates per summand (for the forced-zero rule).
    pub summand_coords: Vec<Vec<usize>>,
    /// Weight cap per summand, if any.
    pub weight_cap: Option<u32>,
}

#[derive(Default, Clone, Debug)]
pub struct RunStats {
    pub calls: u64,
    pub tame: u64,
    pub wild: u64,
    pub budget: u64,
    pub inconsistent: u64,
    pub max_depth: usize,
    pub tame_depths: Vec<usize>,
    pub xor_total: u64,
    pub rows_max: usize,
    pub cols_max: usize,
    pub max_degree: u32,
    pub basis_max: usize,
    pub pruned: u64,
    pub forced: u64,
    pub exhausted: bool,
    pub unresolved_leaves: u64,
    pub solutions: Vec<[u64; W]>,
    pub candidates_rejected: u64,
    pub matrix_cap_hits: u64,
    pub leaf_direct: u64,
    pub prep_ns: u64,
    pub elim_ns: u64,
    pub add_ns: u64,
    pub rref_ns: u64,
    pub subst_ns: u64,
    pub post_ns: u64,
    pub branch_ns: u64,
}

pub fn multisolve(inst: &Instance, cfg: &Config, verify: &mut dyn FnMut(&[u64; W]) -> bool) -> RunStats {
    let mut st = RunStats::default();
    let ones = vec![0u32; inst.summand_coords.len()];
    let assigned = [0u64; W];
    node(inst, cfg, &inst.eqs, 0, ones, assigned, &mut st, verify);
    st
}

#[allow(clippy::too_many_arguments)]
fn node(
    inst: &Instance,
    cfg: &Config,
    eqs: &[Poly],
    depth: usize,
    ones: Vec<u32>,
    assigned: [u64; W],
    st: &mut RunStats,
    verify: &mut dyn FnMut(&[u64; W]) -> bool,
) -> bool {
    if st.exhausted {
        return false;
    }
    if st.calls >= cfg.max_calls {
        st.exhausted = true;
        return false;
    }
    st.calls += 1;
    st.max_depth = st.max_depth.max(depth);
    let (outcome, f4) = truncated_groebner(eqs, cfg.d_max, cfg.budget_xor, cfg.matrix_cap_words);
    st.matrix_cap_hits += f4.matrix_cap_hits;
    st.prep_ns += f4.prep_ns;
    st.elim_ns += f4.elim_ns;
    st.add_ns += f4.add_ns;
    st.rref_ns += f4.rref_ns;
    st.subst_ns += f4.subst_ns;
    st.post_ns += f4.post_ns;
    st.xor_total += f4.xor_words;
    st.rows_max = st.rows_max.max(f4.rows_max);
    st.cols_max = st.cols_max.max(f4.cols_max);
    st.max_degree = st.max_degree.max(f4.max_degree_built);
    st.basis_max = st.basis_max.max(f4.basis_len);
    let wild = match outcome {
        Outcome::Inconsistent => {
            st.tame += 1;
            st.inconsistent += 1;
            st.tame_depths.push(depth);
            return false;
        }
        Outcome::Linear(rref) => match linear_solutions(&rref, inst.n_vars, cfg.max_free) {
            Some(sols) => {
                st.tame += 1;
                st.tame_depths.push(depth);
                for s in sols {
                    let mut pt = s;
                    for k in 0..W {
                        pt[k] |= assigned[k];
                    }
                    if inst.eqs.iter().all(|e| !e.eval(&pt)) && verify(&pt) {
                        st.solutions.push(pt);
                        return true;
                    }
                    st.candidates_rejected += 1;
                }
                return false;
            }
            None => true,
        },
        Outcome::Wild => true,
        Outcome::Budget => {
            st.budget += 1;
            true
        }
    };
    debug_assert!(wild);
    st.wild += 1;
    if depth >= inst.branch.len() {
        // Every coordinate is assigned: the candidate is decided directly,
        // whatever the algebra managed within its budget.
        st.unresolved_leaves += 1;
        if verify(&assigned) {
            st.leaf_direct += 1;
            st.solutions.push(assigned);
            return true;
        }
        return false;
    }
    let (v, s) = inst.branch[depth];
    for b in [false, true] {
        if b {
            if let Some(w) = inst.weight_cap {
                if ones[s] >= w {
                    st.pruned += 1;
                    continue;
                }
            }
        }
        let t_b = std::time::Instant::now();
        let mut next: Vec<Poly> = eqs.iter().map(|e| e.substitute(v, b)).filter(|e| !e.is_zero()).collect();
        st.branch_ns += t_b.elapsed().as_nanos() as u64;
        let mut ones2 = ones.clone();
        let mut assigned2 = assigned;
        let mut depth2 = depth + 1;
        if b {
            ones2[s] += 1;
            assigned2[v / 64] |= 1 << (v % 64);
            if let Some(w) = inst.weight_cap {
                if ones2[s] == w {
                    // force the remaining coordinates of this summand to 0
                    while depth2 < inst.branch.len() && inst.branch[depth2].1 == s {
                        let (v2, _) = inst.branch[depth2];
                        next = next.iter().map(|e| e.substitute(v2, false)).filter(|e| !e.is_zero()).collect();
                        depth2 += 1;
                        st.forced += 1;
                    }
                }
            }
        }
        if node(inst, cfg, &next, depth2, ones2, assigned2, st, verify) {
            return true;
        }
        if st.exhausted {
            return false;
        }
    }
    false
}
