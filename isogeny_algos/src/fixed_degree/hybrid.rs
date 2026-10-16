//! The hybrid approach of the slides ("More guessing"): for smooth d = l^e whose epsilon is too
//! large for the norm-equation solvers alone,
//!
//! 1. guess an l^{e_1}-isogeny phi_1 : E_1 -> E (here: a non-backtracking chain of l-steps on the
//!    quaternion side, `endring::sub_ideals`, which is exactly "use phi_1 to compute End(E) from
//!    End(E_1)": the sub-ideal K of I_1 has O_R(K) = End(E));
//! 2. solve the fixed-degree problem between E and E_2 for degree l^{e - e_1}
//!    (`pipeline::solve_degree` on a Kirschmer–Voight connecting ideal), or guess again;
//! 3. compose: J = J_1 J_2 with J_1 = conj(I_1) K / N(I_1) the ideal of phi_1 (left order
//!    O_1, norm l^{e_1}) and J_2 the solution (left order O_R(K), norm l^{e - e_1}).
//!
//! The result is a left O_1-ideal of norm l^e whose Deuring image is E_2, to be turned into an
//! isogeny by `represent::ideal_to_isogeny`.
use super::connecting::connecting_ideal;
use super::endring::sub_ideals;
use super::pipeline::{solve_degree, Solver};
use crate::int::Int;
use crate::quat::{Alg, Lattice};

#[derive(Clone, Debug)]
pub struct HybridResult {
    /// left O_R(I_1)-ideal of norm l^e
    pub ideal: Lattice,
    /// guessed prefixes tried (each a different E)
    pub guesses: usize,
    pub e1: u32,
}

/// All non-backtracking chains of `e1` l-steps from the left O_0-ideal `i1` (as the final
/// sub-ideals; distinct chains give distinct ideals).
pub fn prefixes(alg: &Alg, o0: &Lattice, i1: &Lattice, ell: u64, e1: u32) -> Vec<Lattice> {
    let mut layer: Vec<(Lattice, Option<Lattice>)> = vec![(i1.clone(), None)];
    for _ in 0..e1 {
        let mut next = vec![];
        for (k, prev) in &layer {
            let back = prev.as_ref().map(|pv| pv.scale(&Int::from(ell), &Int::one()));
            for c in sub_ideals(alg, o0, k, ell) {
                if Some(&c) == back.as_ref() {
                    continue;
                }
                next.push((c, Some(k.clone())));
            }
        }
        layer = next;
    }
    layer.into_iter().map(|(k, _)| k).collect()
}

/// Hybrid search for an ideal of norm l^e from O_1 = O_R(i1) to O_2 = O_R(i2).
pub fn hybrid_search(
    alg: &Alg,
    o0: &Lattice,
    i1: &Lattice,
    i2: &Lattice,
    ell: u64,
    e: u32,
    e1: u32,
    solver: Solver,
) -> Option<HybridResult> {
    let o1 = i1.right_order(alg, o0);
    let o2 = i2.right_order(alg, o0);
    let n1 = i1.norm_int(o0);
    let rest = Int::from(ell).pow(e - e1);
    let mut guesses = 0usize;
    for k in prefixes(alg, o0, i1, ell, e1) {
        guesses += 1;
        let ok = k.right_order(alg, o0);
        let conn = connecting_ideal(alg, &ok, &o2);
        let Some(sol) = solve_degree(alg, &ok, &conn, &rest, solver, None) else {
            continue;
        };
        let j1 = i1.conj().mul(alg, &k).scale(&Int::one(), &n1);
        assert!(j1.norm(&o1) == (Int::from(ell).pow(e1), Int::one()));
        let jd = j1.mul(alg, &sol.ideal);
        assert!(jd.norm(&o1) == (Int::from(ell).pow(e), Int::one()));
        return Some(HybridResult {
            ideal: jd,
            guesses,
            e1,
        });
    }
    None
}
