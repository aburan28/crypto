//! **Split and hybrid point-decomposition systems**: the solving degree of
//! the summation-polynomial PDP when one decomposition `R = P₁ + … + P_m`
//! is cut into blocks joined by full-field auxiliary abscissae, and when
//! `g` Boolean unknowns are guessed (Bettale–Faugère–Perret hybrid).
//!
//! Study: `research/pdp_split_dreg_20261011/`.  A plan is a list of block
//! sizes summing to `m`.  Block `i` is one summation-polynomial link:
//!
//! - one block `[m]` is the unsplit `S_{m+1}(x₁, …, x_m, x_R) = 0`;
//! - the first of several blocks is `S_{s₁+1}(x…, u₁) = 0`;
//! - a middle block is `S_{sᵢ+2}(u_{i−1}, x…, uᵢ) = 0`;
//! - the last block is `S_{s_k+2}(u_{k−1}, x…, x_R) = 0`.
//!
//! Each `uⱼ` is an unknown abscissa over all of `F_{2^n}` (`n` Boolean
//! unknowns after Weil descent); each `xᵢ` ranges over the factor-base
//! subspace `V` (`ℓ` unknowns).  Links are built from the repository's
//! symbolic `S₃` ([`sym_semaev_s3`], any `b`) and `S₄`
//! ([`sym_semaev_s4`], `b = 1`, last argument a known constant), so the
//! supported plans are `[3]`, `[2, 1]`, `[2, 2]` and `[2, 1, …, 1]`; the
//! unsplit `m = 4` system needs a symbolic binary `S₅`, which the
//! repository does not have, and is refused rather than approximated.
//!
//! Every solve is checked twice, independently of the algebra: each
//! Boolean solution's summands are lifted to curve points and some signed
//! sum must have abscissa `x_R`, and the verified set is compared with an
//! exhaustive enumeration of the decompositions of `R` over `V`.

use std::collections::BTreeSet;
use std::time::{Duration, Instant};

use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::ic_boundary::BinaryInstance;
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_groebner::{
    sym_semaev_s3, sym_semaev_s4, FieldStructure, SymElement, MAX_VARS,
};
use crate::cryptanalysis::pq_f4_f2::{groebner_basis_f4, solutions_from_reduced_basis, F4Stats};
use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

/// A Weil-descended decomposition system for one target abscissa.
#[derive(Clone, Debug)]
pub struct SplitSystem {
    pub plan: Vec<usize>,
    pub m: usize,
    pub ell: usize,
    pub n: u32,
    /// `m·ℓ` summand unknowns, then `(blocks − 1)·n` auxiliary ones.
    pub n_vars: usize,
    pub equations: Vec<F2BoolPoly>,
}

impl SplitSystem {
    /// Highest total degree of any input equation.
    pub fn input_degree(&self) -> u32 {
        self.equations
            .iter()
            .flat_map(|p| p.terms.iter().map(|t| t.degree()))
            .max()
            .unwrap_or(0)
    }

    /// Terms over all input equations.
    pub fn input_terms(&self) -> usize {
        self.equations.iter().map(|p| p.terms.len()).sum()
    }
}

/// The plan's short name: `unsplit` for one block, else `split-2+1+1`.
pub fn plan_label(plan: &[usize]) -> String {
    if plan.len() == 1 {
        "unsplit".into()
    } else {
        let parts: Vec<String> = plan.iter().map(|s| s.to_string()).collect();
        format!("split-{}", parts.join("+"))
    }
}

/// The standard factor-base subspace `span(1, z, …, z^{ℓ−1})`, the base
/// of ecbench's `koblitz-standard-subspace:dimension=ℓ`.
pub fn standard_basis(ell: usize, n: u32) -> Vec<F2mElement> {
    (0..ell as u32)
        .map(|t| F2mElement::from_bit_positions(&[t], n))
        .collect()
}

/// Build the system of `plan` for target abscissa `x_r` over the
/// subspace spanned by `basis`, on a curve with constant `b`.
pub fn build_split_system(
    plan: &[usize],
    basis: &[F2mElement],
    x_r: &F2mElement,
    b: &F2mElement,
    st: &FieldStructure,
) -> Result<SplitSystem, String> {
    let n = st.n;
    let ell = basis.len();
    let m: usize = plan.iter().sum();
    if plan.is_empty() || plan.contains(&0) {
        return Err("plan needs nonzero block sizes".into());
    }
    let k = plan.len();
    let n_vars = m * ell + (k - 1) * n as usize;
    if n_vars > MAX_VARS {
        return Err(format!(
            "{n_vars} Boolean unknowns exceed the cap {MAX_VARS}"
        ));
    }
    let b_is_one = *b == F2mElement::one(n);
    let summands: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(basis, i * ell, n, n_vars))
        .collect();
    let aux: Vec<SymElement> = (0..k - 1)
        .map(|j| SymElement::from_free_vars(m * ell + j * n as usize, n, n_vars))
        .collect();
    let target = SymElement::constant(x_r, n, n_vars);

    let mut equations = Vec::new();
    let mut next = 0usize;
    for (i, &size) in plan.iter().enumerate() {
        let xs = &summands[next..next + size];
        next += size;
        let first = i == 0;
        let last = i + 1 == k;
        // The link's symbolic arguments, and whether it closes on x_R.
        let mut args: Vec<&SymElement> = Vec::new();
        if !first {
            args.push(&aux[i - 1]);
        }
        args.extend(xs.iter());
        match (args.len(), last) {
            // S₃(a, c, d) with d = u_i or x_R.
            (2, _) => {
                let close = if last { &target } else { &aux[i] };
                equations.extend(sym_semaev_s3(args[0], args[1], close, b, st));
            }
            // S₄(a, c, d, x_R): only a known constant can close it.
            (3, true) => {
                if !b_is_one {
                    return Err("the symbolic S4 link is specialised to b = 1".into());
                }
                equations.extend(sym_semaev_s4(args[0], args[1], args[2], x_r, st));
            }
            (l, _) => {
                return Err(format!(
                    "link {i} needs S_{} with a free last argument; plan {plan:?} is not \
                     supported (S3 links anywhere, S4 only as the closing link)",
                    l + 1
                ));
            }
        }
    }
    equations.retain(|p| !p.is_zero());
    Ok(SplitSystem {
        plan: plan.to_vec(),
        m,
        ell,
        n,
        n_vars,
        equations,
    })
}

/// One hybrid branch: the first `g` unknowns fixed to the bits of `guess`.
#[derive(Clone, Debug, serde::Serialize)]
pub struct BranchOutcome {
    pub guess: u64,
    pub timed_out: bool,
    /// Boolean solutions (full assignments), when the basis resolved them.
    pub solutions: Option<Vec<u64>>,
    pub solving_degree: u32,
    pub degree_reached: u32,
    pub max_poly_degree: u32,
    pub matrix_rows_max: u64,
    pub matrix_cols_max: u64,
    pub word_xors: u64,
    pub peak_matrix_bytes: u64,
    pub steps: u64,
    pub wall_ns: u64,
}

/// Specialise the first `g` unknowns to `guess` and adjoin their linear
/// pins, so the reduced basis still determines every unknown.
pub fn specialise(sys: &SplitSystem, g: usize, guess: u64) -> Vec<F2BoolPoly> {
    let mut eqs: Vec<F2BoolPoly> = sys
        .equations
        .iter()
        .map(|p| {
            let mut q = p.clone();
            for v in 0..g {
                q = q.substitute(v as u32, (guess >> v) & 1 == 1);
            }
            q
        })
        .filter(|p| !p.is_zero())
        .collect();
    for v in 0..g {
        let mut monos = vec![F2BoolMono::var(v as u32)];
        if (guess >> v) & 1 == 1 {
            monos.push(F2BoolMono::one());
        }
        eqs.push(F2BoolPoly::from_monos(monos, sys.n_vars));
    }
    eqs
}

/// Solve one branch with the repository's F4 over `F₂[v]/(v² − v)`.
pub fn solve_branch(
    sys: &SplitSystem,
    g: usize,
    guess: u64,
    budget: Option<Duration>,
    max_free: usize,
) -> BranchOutcome {
    let eqs = specialise(sys, g, guess);
    let (gb, st): (Vec<F2BoolPoly>, F4Stats) = groebner_basis_f4(eqs, sys.n_vars, budget);
    let solutions = if st.timed_out {
        None
    } else {
        let mut tests = 0u64;
        solutions_from_reduced_basis(&gb, sys.n_vars, max_free, &mut tests)
    };
    BranchOutcome {
        guess,
        timed_out: st.timed_out,
        solutions,
        solving_degree: st.solving_degree,
        degree_reached: st.degree_reached,
        max_poly_degree: st.max_poly_degree,
        matrix_rows_max: st.matrix_rows_max,
        matrix_cols_max: st.matrix_cols_max,
        word_xors: st.word_xors,
        peak_matrix_bytes: st.peak_matrix_bytes,
        steps: st.steps,
        wall_ns: st.wall_ns,
    }
}

/// A decomposition up to order and sign: the summands' abscissae, sorted.
pub type Decomposition = Vec<u64>;

/// Lift an abscissa to a curve point (`None` if it lies on the twist).
pub fn lift_x(inst: &BinaryInstance, x: u64) -> Option<FastPoint> {
    let gf = &inst.gf;
    if x == 0 {
        // y² = b; squaring is a bijection, so y = b^{2^{n−1}}.
        return Some(FastPoint::affine(0, gf.sqr_k(inst.b, inst.n - 1)));
    }
    let xi = gf.inv(x);
    let c = x ^ inst.a ^ gf.mul(inst.b, gf.sqr(xi));
    let w = inst.artin_schreier.solve(c)?;
    let p = FastPoint::affine(x, gf.mul(x, w));
    debug_assert!(inst.fast.is_on_curve(p));
    Some(p)
}

/// Does some signed sum of the lifted `xs` have abscissa `x_r`?  `None`
/// when a summand does not lift.
pub fn verify_decomposition(inst: &BinaryInstance, xs: &[u64], x_r: u64) -> Option<bool> {
    let pts: Vec<FastPoint> = xs.iter().map(|&x| lift_x(inst, x)).collect::<Option<_>>()?;
    // The first sign is free: ±R are both accepted.
    let k = pts.len();
    for signs in 0..(1u64 << (k - 1)) {
        let mut acc = pts[0];
        for (i, p) in pts.iter().enumerate().skip(1) {
            let q = if (signs >> (i - 1)) & 1 == 1 {
                inst.fast.neg(*p)
            } else {
                *p
            };
            acc = inst.fast.add(acc, q);
        }
        if !acc.infinity && acc.x == x_r {
            return Some(true);
        }
    }
    Some(false)
}

/// The summands' abscissae of a Boolean assignment, in the field's word
/// representation (the standard basis puts bit `t` of block `i` at `z^t`).
pub fn summands_of(sys: &SplitSystem, basis_words: &[u64], assignment: u64) -> Vec<u64> {
    (0..sys.m)
        .map(|i| {
            let mut x = 0u64;
            for (t, w) in basis_words.iter().enumerate() {
                if (assignment >> (i * sys.ell + t)) & 1 == 1 {
                    x ^= w;
                }
            }
            x
        })
        .collect()
}

/// Every decomposition of `R` over `V` with `m` summands, by exhaustive
/// enumeration of multisets of liftable abscissae — the ground truth the
/// algebraic solves are checked against.
pub fn exhaustive_decompositions(
    inst: &BinaryInstance,
    basis_words: &[u64],
    m: usize,
    x_r: u64,
) -> BTreeSet<Decomposition> {
    let ell = basis_words.len();
    let liftable: Vec<u64> = (0..1u64 << ell)
        .map(|c| {
            basis_words
                .iter()
                .enumerate()
                .filter(|(t, _)| (c >> t) & 1 == 1)
                .fold(0u64, |a, (_, w)| a ^ w)
        })
        .filter(|&x| lift_x(inst, x).is_some())
        .collect();
    let mut out = BTreeSet::new();
    let mut idx = vec![0usize; m];
    loop {
        let xs: Vec<u64> = idx.iter().map(|&i| liftable[i]).collect();
        if verify_decomposition(inst, &xs, x_r) == Some(true) {
            let mut d = xs;
            d.sort_unstable();
            out.insert(d);
        }
        // Next non-decreasing index tuple.
        let mut p = m;
        loop {
            if p == 0 {
                return out;
            }
            p -= 1;
            if idx[p] + 1 < liftable.len() {
                idx[p] += 1;
                for q in p + 1..m {
                    idx[q] = idx[p];
                }
                break;
            }
        }
    }
}

/// Liftable abscissae in `V` (points up to sign): the factor-base size.
pub fn factor_base_size(inst: &BinaryInstance, basis_words: &[u64]) -> usize {
    let ell = basis_words.len();
    (0..1u64 << ell)
        .filter(|&c| {
            let x = basis_words
                .iter()
                .enumerate()
                .filter(|(t, _)| (c >> t) & 1 == 1)
                .fold(0u64, |a, (_, w)| a ^ w);
            lift_x(inst, x).is_some()
        })
        .count()
}

/// Everything one (plan, guess-count) solve of one target recorded.
#[derive(Clone, Debug, serde::Serialize)]
pub struct HybridRun {
    pub plan: String,
    pub g: usize,
    pub n_vars: usize,
    pub n_equations: usize,
    pub input_degree: u32,
    pub input_terms: usize,
    pub build_ns: u64,
    pub branches: usize,
    pub branches_timed_out: usize,
    pub branches_unresolved: usize,
    /// Max over branches, and the mean.
    pub solving_degree_max: u32,
    pub solving_degree_mean: f64,
    pub degree_reached_max: u32,
    pub matrix_rows_max: u64,
    pub matrix_cols_max: u64,
    pub peak_matrix_bytes: u64,
    /// Sums over branches.
    pub word_xors: u64,
    pub f4_wall_ns: u64,
    /// Boolean solutions found, and how they checked out on the curve.
    pub algebraic_solutions: usize,
    pub unliftable_solutions: usize,
    pub false_solutions: usize,
    pub verified: Vec<Decomposition>,
    /// The verified set equals the exhaustive one (`None` if incomplete).
    pub complete_and_exact: Option<bool>,
    pub per_branch: Vec<BranchOutcome>,
}

/// Solve `sys` with `g` guessed unknowns, enumerating every branch.
#[allow(clippy::too_many_arguments)]
pub fn run_hybrid(
    inst: &BinaryInstance,
    sys: &SplitSystem,
    basis_words: &[u64],
    x_r: u64,
    g: usize,
    budget: Option<Duration>,
    max_free: usize,
    truth: &BTreeSet<Decomposition>,
    build_ns: u64,
) -> HybridRun {
    let mut per_branch = Vec::with_capacity(1 << g);
    let mut verified = BTreeSet::new();
    let (mut algebraic, mut unliftable, mut false_sols) = (0usize, 0usize, 0usize);
    for guess in 0..(1u64 << g) {
        let bo = solve_branch(sys, g, guess, budget, max_free);
        if let Some(sols) = &bo.solutions {
            for &a in sols {
                algebraic += 1;
                let xs = summands_of(sys, basis_words, a);
                match verify_decomposition(inst, &xs, x_r) {
                    None => unliftable += 1,
                    Some(false) => false_sols += 1,
                    Some(true) => {
                        let mut d = xs;
                        d.sort_unstable();
                        verified.insert(d);
                    }
                }
            }
        }
        per_branch.push(bo);
    }
    let timed_out = per_branch.iter().filter(|b| b.timed_out).count();
    let unresolved = per_branch
        .iter()
        .filter(|b| !b.timed_out && b.solutions.is_none())
        .count();
    let complete = timed_out == 0 && unresolved == 0;
    let sd: Vec<u32> = per_branch.iter().map(|b| b.solving_degree).collect();
    HybridRun {
        plan: plan_label(&sys.plan),
        g,
        n_vars: sys.n_vars,
        n_equations: sys.equations.len(),
        input_degree: sys.input_degree(),
        input_terms: sys.input_terms(),
        build_ns,
        branches: per_branch.len(),
        branches_timed_out: timed_out,
        branches_unresolved: unresolved,
        solving_degree_max: sd.iter().copied().max().unwrap_or(0),
        solving_degree_mean: sd.iter().map(|&d| d as f64).sum::<f64>() / sd.len() as f64,
        degree_reached_max: per_branch
            .iter()
            .map(|b| b.degree_reached)
            .max()
            .unwrap_or(0),
        matrix_rows_max: per_branch
            .iter()
            .map(|b| b.matrix_rows_max)
            .max()
            .unwrap_or(0),
        matrix_cols_max: per_branch
            .iter()
            .map(|b| b.matrix_cols_max)
            .max()
            .unwrap_or(0),
        peak_matrix_bytes: per_branch
            .iter()
            .map(|b| b.peak_matrix_bytes)
            .max()
            .unwrap_or(0),
        word_xors: per_branch.iter().map(|b| b.word_xors).sum(),
        f4_wall_ns: per_branch.iter().map(|b| b.wall_ns).sum(),
        algebraic_solutions: algebraic,
        unliftable_solutions: unliftable,
        false_solutions: false_sols,
        complete_and_exact: complete.then(|| &verified == truth),
        verified: verified.into_iter().collect(),
        per_branch,
    }
}

/// Build `plan`'s system for `x_r` on `inst` over the standard subspace of
/// dimension `ell`, timing the construction.
pub fn build_for_instance(
    inst: &BinaryInstance,
    plan: &[usize],
    ell: usize,
    x_r: u64,
) -> Result<(SplitSystem, u64), String> {
    let n = inst.n;
    let st = FieldStructure::new(n, &inst.irreducible);
    let basis = standard_basis(ell, n);
    let b = inst.gf.to_element(inst.b);
    let xr = inst.gf.to_element(x_r);
    let t0 = Instant::now();
    let sys = build_split_system(plan, &basis, &xr, &b, &st)?;
    Ok((sys, t0.elapsed().as_nanos() as u64))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::koblitz_instance;

    /// Unsplit and split systems, with and without guessing, find exactly
    /// the decompositions exhaustive enumeration finds.
    #[test]
    fn plans_agree_with_exhaustive_m3() {
        let inst = koblitz_instance(1, 11).expect("K_1 / GF(2^11)");
        let ell = 4;
        let words: Vec<u64> = (0..ell).map(|t| 1u64 << t).collect();
        let mut hits = 0;
        for k in 1..40u64 {
            let r = inst.fast.mul_u64(inst.generator, k * 7919 + 3);
            if r.infinity {
                continue;
            }
            let truth = exhaustive_decompositions(&inst, &words, 3, r.x);
            hits += usize::from(!truth.is_empty());
            for plan in [vec![3], vec![2, 1]] {
                for g in [0usize, 2] {
                    let (sys, ns) = build_for_instance(&inst, &plan, ell, r.x).unwrap();
                    let run = run_hybrid(&inst, &sys, &words, r.x, g, None, 20, &truth, ns);
                    assert_eq!(run.false_solutions, 0);
                    assert_eq!(run.complete_and_exact, Some(true), "{plan:?} g={g} k={k}");
                }
            }
            if hits >= 3 {
                break;
            }
        }
        assert!(hits > 0, "no decomposable target among the samples");
    }

    /// About two minutes in release (the `[2, 2]` link is an `S₄` with a
    /// free argument): run with `--ignored`.
    #[test]
    #[ignore = "slow: about two minutes in release"]
    fn plans_agree_with_exhaustive_m4() {
        let inst = koblitz_instance(1, 11).expect("K_1 / GF(2^11)");
        let ell = 3;
        let words: Vec<u64> = (0..ell).map(|t| 1u64 << t).collect();
        let r = inst.fast.mul_u64(inst.generator, 12345);
        let truth = exhaustive_decompositions(&inst, &words, 4, r.x);
        for plan in [vec![2, 2], vec![2, 1, 1]] {
            let (sys, ns) = build_for_instance(&inst, &plan, ell, r.x).unwrap();
            let run = run_hybrid(&inst, &sys, &words, r.x, 0, None, 20, &truth, ns);
            assert_eq!(run.false_solutions, 0);
            assert_eq!(run.complete_and_exact, Some(true), "{plan:?}");
        }
    }

    #[test]
    fn unsupported_plans_are_refused() {
        let inst = koblitz_instance(1, 11).unwrap();
        assert!(build_for_instance(&inst, &[4], 3, 5).is_err());
        assert!(build_for_instance(&inst, &[3, 1], 3, 5).is_err());
    }
}
