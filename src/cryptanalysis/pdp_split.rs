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
use crate::cryptanalysis::pq_f4_f2::{groebner_basis_f4, F4Stats};
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
    let m: usize = plan.iter().sum();
    let bases = vec![basis.to_vec(); m];
    build_split_system_bases(plan, &bases, x_r, b, st)
}

/// [`build_split_system`] with summand `i` ranging over the span of
/// `bases[i]` — the τ-slot form: over a Frobenius closure, summand `i`
/// lives in the rotation `σ^{jᵢ}(V)`.  Every basis has the same length `ℓ`.
pub fn build_split_system_bases(
    plan: &[usize],
    bases: &[Vec<F2mElement>],
    x_r: &F2mElement,
    b: &F2mElement,
    st: &FieldStructure,
) -> Result<SplitSystem, String> {
    let n = st.n;
    let m: usize = plan.iter().sum();
    if plan.is_empty() || plan.contains(&0) {
        return Err("plan needs nonzero block sizes".into());
    }
    if bases.len() != m || bases.iter().any(|b| b.len() != bases[0].len()) {
        return Err(format!("need {m} bases of one dimension"));
    }
    let ell = bases[0].len();
    let k = plan.len();
    let n_vars = m * ell + (k - 1) * n as usize;
    // `bitop/v1`: the symbolic build is charged by `SymElement` and the
    // Semaev links; here, the scan that drops zero equations
    crate::cryptanalysis::bitop_tally::word_ops(n as u64 * plan.len() as u64);
    if n_vars > MAX_VARS {
        return Err(format!(
            "{n_vars} Boolean unknowns exceed the cap {MAX_VARS}"
        ));
    }
    let b_is_one = *b == F2mElement::one(n);
    let summands: Vec<SymElement> = (0..m)
        .map(|i| SymElement::from_subspace_vars(&bases[i], i * ell, n, n_vars))
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
    /// The XORs the eliminations performed, tables included (`F4Stats`).
    pub word_xors_performed: u64,
    pub peak_matrix_bytes: u64,
    pub steps: u64,
    pub wall_ns: u64,
}

/// Specialise the first `g` unknowns to `guess` and adjoin their linear
/// pins, so the reduced basis still determines every unknown.
pub fn specialise(sys: &SplitSystem, g: usize, guess: u64) -> Vec<F2BoolPoly> {
    // `bitop/v1`: `g` substitutions of every equation (each at most its
    // terms), and the pins
    crate::cryptanalysis::bitop_tally::word_ops(
        sys.equations
            .iter()
            .map(|p| {
                1 + g as u64 * crate::cryptanalysis::bitop_tally::poly::substitute(p.terms.len())
            })
            .sum::<u64>()
            + 3 * g as u64,
    );
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

/// The `F₂` points of the ideal `basis` generates, for any generating set.
///
/// `pq_f4_f2::solutions_from_reduced_basis` assumes every linear element's
/// tail holds free variables only.  F4's output does not always satisfy
/// that (study amendment 1: a basis whose linear elements reference each
/// other's leading variables read as no solutions).  Here the linear
/// elements are brought to reduced row-echelon form first, every
/// assignment of the remaining free variables is completed through them,
/// and a point is kept only if every element of `basis` vanishes on it.
/// `None` above `max_free` free variables.
pub fn solutions_from_basis(
    basis: &[F2BoolPoly],
    n_vars: usize,
    max_free: usize,
) -> Option<Vec<u64>> {
    // `bitop/v1`: the linearity test and fold of every term, one word
    // operation a pivot test of the echelon, and per enumerated point one
    // a free variable and pivot, and one a term evaluated
    let mut ops: u64 = basis.iter().map(|p| 2 * p.terms.len() as u64).sum();
    let mut tally = crate::cryptanalysis::bitop_tally::Pending::default();
    // A linear row: variable mask and constant bit.
    let mut rows: Vec<(u64, bool)> = basis
        .iter()
        .filter(|p| p.terms.iter().all(|t| t.degree() <= 1))
        .map(|p| {
            let mut mask = 0u64;
            let mut c = false;
            for t in &p.terms {
                if t.mask == 0 {
                    c = !c;
                } else {
                    mask ^= t.mask;
                }
            }
            (mask, c)
        })
        .collect();
    // Reduced row-echelon form, pivoting on the highest variable of a row.
    let mut pivots: Vec<(u32, u64, bool)> = Vec::new();
    for i in 0..rows.len() {
        let (mut mask, mut c) = rows[i];
        ops += 2 * pivots.len() as u64 + 1;
        for &(pv, pm, pc) in &pivots {
            if (mask >> pv) & 1 == 1 {
                mask ^= pm;
                c ^= pc;
            }
        }
        if mask == 0 {
            if c {
                return Some(Vec::new()); // 1 = 0: no points
            }
            continue;
        }
        let pv = 63 - mask.leading_zeros();
        for p in pivots.iter_mut() {
            if (p.1 >> pv) & 1 == 1 {
                p.1 ^= mask;
                p.2 ^= c;
            }
        }
        pivots.push((pv, mask, c));
        rows[i] = (mask, c);
    }
    tally.word_ops += ops;
    let pivot_mask: u64 = pivots.iter().fold(0, |a, p| a | (1u64 << p.0));
    let free: Vec<usize> = (0..n_vars)
        .filter(|&v| (pivot_mask >> v) & 1 == 0)
        .collect();
    if free.len() > max_free {
        return None;
    }
    let mut out = Vec::new();
    for a in 0u64..(1u64 << free.len()) {
        let mut point = 0u64;
        for (bit, &v) in free.iter().enumerate() {
            if (a >> bit) & 1 == 1 {
                point |= 1u64 << v;
            }
        }
        // Each pivot row is v_p + Σ(free) + c with no other pivot in it.
        for &(pv, pm, pc) in &pivots {
            let tail = pm & !(1u64 << pv);
            let val = ((point & tail).count_ones() & 1 == 1) ^ pc;
            if val {
                point |= 1u64 << pv;
            }
        }
        tally.word_ops += (free.len() + pivots.len()) as u64;
        if basis.iter().all(|g| {
            tally.word_ops += g.terms.len() as u64;
            g.eval(point) == 0
        }) {
            out.push(point);
        }
    }
    Some(out)
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
        solutions_from_basis(&gb, sys.n_vars, max_free)
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
        word_xors_performed: st.word_xors_performed,
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
    summands_of_bases(sys, &vec![basis_words.to_vec(); sys.m], assignment)
}

/// [`summands_of`] with summand `i` over the span of `bases_words[i]`.
pub fn summands_of_bases(sys: &SplitSystem, bases_words: &[Vec<u64>], assignment: u64) -> Vec<u64> {
    (0..sys.m)
        .map(|i| {
            let mut x = 0u64;
            for (t, w) in bases_words[i].iter().enumerate() {
                if (assignment >> (i * sys.ell + t)) & 1 == 1 {
                    x ^= w;
                }
            }
            x
        })
        .collect()
}

/// Every element of the `F₂`-span of `basis_words`.
pub fn span_words(basis_words: &[u64]) -> Vec<u64> {
    (0..1u64 << basis_words.len())
        .map(|c| {
            basis_words
                .iter()
                .enumerate()
                .filter(|(t, _)| (c >> t) & 1 == 1)
                .fold(0u64, |a, (_, w)| a ^ w)
        })
        .collect()
}

/// The basis of the rotation `σ^j(V)`: each basis word raised to `2^j`.
pub fn rotate_words(inst: &BinaryInstance, basis_words: &[u64], j: u32) -> Vec<u64> {
    basis_words.iter().map(|&w| inst.gf.sqr_k(w, j)).collect()
}

/// The slot tuples of one target over a Frobenius closure: summand `i`
/// over `σ^{jᵢ}(V)`, with `j₁ = 0` and `j₂ ≤ … ≤ j_m` (the systems are
/// symmetric in their summands, so a multiset of rotations is one system).
/// `[0, 0, …, 0]` is the plain `V^m` system.
pub fn slot_tuples(n: u32, m: usize) -> Vec<Vec<u32>> {
    let mut out = Vec::new();
    let mut cur = vec![0u32; m];
    fn rec(n: u32, m: usize, pos: usize, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if pos == m {
            out.push(cur.clone());
            return;
        }
        let lo = if pos <= 1 { 0 } else { cur[pos - 1] };
        for j in lo..n {
            cur[pos] = j;
            rec(n, m, pos + 1, cur, out);
        }
    }
    if m == 1 {
        return vec![vec![0]];
    }
    rec(n, m, 1, &mut cur, &mut out);
    out
}

/// The Frobenius closure `∪_j σ^j(V)` as a factor base: how many liftable
/// abscissae it holds, and how many columns it needs once each signed
/// Frobenius orbit `{±τ^k P}` is folded to one (an abscissa orbit under
/// squaring is exactly such a point orbit).
#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct ClosureCensus {
    pub abscissae: usize,
    pub columns: usize,
}

pub fn closure_census(inst: &BinaryInstance, basis_words: &[u64]) -> ClosureCensus {
    let mut xs = BTreeSet::new();
    for j in 0..inst.n {
        for x in span_words(&rotate_words(inst, basis_words, j)) {
            if lift_x(inst, x).is_some() {
                xs.insert(x);
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut columns = 0usize;
    for &x in &xs {
        if seen.contains(&x) {
            continue;
        }
        columns += 1;
        let mut y = x;
        for _ in 0..inst.n {
            seen.insert(y);
            y = inst.gf.sqr(y);
        }
    }
    ClosureCensus {
        abscissae: xs.len(),
        columns,
    }
}

/// [`exhaustive_decompositions`] with summand `i` over the span of
/// `bases_words[i]`: every tuple of liftable abscissae, one from each
/// space, whose points sum to `±R` with some signs, as sorted multisets.
pub fn exhaustive_decompositions_bases(
    inst: &BinaryInstance,
    bases_words: &[Vec<u64>],
    x_r: u64,
) -> BTreeSet<Decomposition> {
    let lists: Vec<Vec<u64>> = bases_words
        .iter()
        .map(|b| {
            span_words(b)
                .into_iter()
                .filter(|&x| lift_x(inst, x).is_some())
                .collect()
        })
        .collect();
    let mut out = BTreeSet::new();
    if lists.iter().any(|l| l.is_empty()) {
        return out;
    }
    let m = lists.len();
    let mut idx = vec![0usize; m];
    loop {
        let xs: Vec<u64> = idx.iter().zip(&lists).map(|(&i, l)| l[i]).collect();
        if verify_decomposition(inst, &xs, x_r) == Some(true) {
            let mut d = xs;
            d.sort_unstable();
            out.insert(d);
        }
        let mut p = m;
        loop {
            if p == 0 {
                return out;
            }
            p -= 1;
            if idx[p] + 1 < lists[p].len() {
                idx[p] += 1;
                for q in idx.iter_mut().skip(p + 1) {
                    *q = 0;
                }
                break;
            }
        }
    }
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
    let bases = vec![basis_words.to_vec(); sys.m];
    run_hybrid_bases(inst, sys, &bases, x_r, g, budget, max_free, truth, build_ns)
}

/// [`run_hybrid`] with summand `i` over the span of `bases_words[i]`.
#[allow(clippy::too_many_arguments)]
pub fn run_hybrid_bases(
    inst: &BinaryInstance,
    sys: &SplitSystem,
    bases_words: &[Vec<u64>],
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
                let xs = summands_of_bases(sys, bases_words, a);
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

/// Build `plan`'s system for `x_r` with summand `i` over `σ^{slots[i]}(V)`,
/// `V` the standard subspace of dimension `ell`; returns the system, the
/// per-summand basis words, and the construction time.
pub fn build_for_instance_slots(
    inst: &BinaryInstance,
    plan: &[usize],
    ell: usize,
    slots: &[u32],
    x_r: u64,
) -> Result<(SplitSystem, Vec<Vec<u64>>, u64), String> {
    let n = inst.n;
    let st = FieldStructure::new(n, &inst.irreducible);
    let v: Vec<u64> = (0..ell).map(|t| 1u64 << t).collect();
    let words: Vec<Vec<u64>> = slots.iter().map(|&j| rotate_words(inst, &v, j)).collect();
    let bases: Vec<Vec<F2mElement>> = words
        .iter()
        .map(|w| w.iter().map(|&x| inst.gf.to_element(x)).collect())
        .collect();
    let b = inst.gf.to_element(inst.b);
    let xr = inst.gf.to_element(x_r);
    let t0 = Instant::now();
    let sys = build_split_system_bases(plan, &bases, &xr, &b, &st)?;
    Ok((sys, words, t0.elapsed().as_nanos() as u64))
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

    /// Amendment 1's regression: the S₃ chain at `K_0 / GF(2^13)`, `ℓ = 5`,
    /// target abscissa 4035, `g = 0`.  F4's basis there has linear elements
    /// that reference each other's leading variables; the decomposition
    /// `{5, 21, 22}` must still be read out.
    #[test]
    fn non_interreduced_basis_still_yields_points() {
        let inst = koblitz_instance(0, 13).expect("K_0 / GF(2^13)");
        let words: Vec<u64> = (0..5).map(|t| 1u64 << t).collect();
        let x_r = 4035;
        let truth = exhaustive_decompositions(&inst, &words, 3, x_r);
        assert_eq!(truth.len(), 1);
        let (sys, ns) = build_for_instance(&inst, &[2, 1], 5, x_r).unwrap();
        let run = run_hybrid(&inst, &sys, &words, x_r, 0, None, 24, &truth, ns);
        assert_eq!(run.complete_and_exact, Some(true));
        assert_eq!(run.verified, vec![vec![5, 21, 22]]);
    }

    /// Over rotated per-summand bases, unsplit and split systems (with and
    /// without guessing) find exactly the exhaustive decompositions, and
    /// the all-zero slot reproduces the plain system.
    #[test]
    fn slot_systems_agree_with_exhaustive() {
        let inst = koblitz_instance(1, 11).expect("K_1 / GF(2^11)");
        let ell = 4;
        let v: Vec<u64> = (0..ell).map(|t| 1u64 << t).collect();
        assert_eq!(slot_tuples(11, 3).len(), 66);
        let mut hits = 0;
        for (k, slots) in [[0u32, 0, 0], [0, 1, 5], [0, 3, 3], [0, 2, 9]]
            .iter()
            .enumerate()
        {
            let r = inst
                .fast
                .mul_u64(inst.generator, 4099 * (k as u64 + 1) + 17);
            let words: Vec<Vec<u64>> = slots.iter().map(|&j| rotate_words(&inst, &v, j)).collect();
            let truth = exhaustive_decompositions_bases(&inst, &words, r.x);
            if *slots == [0, 0, 0] {
                assert_eq!(truth, exhaustive_decompositions(&inst, &v, 3, r.x));
            }
            hits += truth.len();
            for plan in [vec![3], vec![2, 1]] {
                for g in [0usize, 2] {
                    let (sys, w, ns) =
                        build_for_instance_slots(&inst, &plan, ell, slots, r.x).unwrap();
                    let run = run_hybrid_bases(&inst, &sys, &w, r.x, g, None, 24, &truth, ns);
                    assert_eq!(run.false_solutions, 0);
                    assert_eq!(
                        run.complete_and_exact,
                        Some(true),
                        "{plan:?} {slots:?} g={g}"
                    );
                }
            }
        }
        let _ = hits;
    }

    #[test]
    fn closure_folds_by_frobenius() {
        let inst = koblitz_instance(1, 11).unwrap();
        let v: Vec<u64> = (0..4).map(|t| 1u64 << t).collect();
        let c = closure_census(&inst, &v);
        // Every column is a whole squaring orbit of liftable abscissae.
        assert!(c.columns <= c.abscissae);
        assert!(c.abscissae <= 11 * 16);
        assert!(c.columns >= factor_base_size(&inst, &v) / 11);
    }

    #[test]
    fn unsupported_plans_are_refused() {
        let inst = koblitz_instance(1, 11).unwrap();
        assert!(build_for_instance(&inst, &[4], 3, 5).is_err());
        assert!(build_for_instance(&inst, &[3, 1], 3, 5).is_err());
    }
}
