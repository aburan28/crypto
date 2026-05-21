#![allow(dead_code)]
// F4 Gröbner basis algorithm — STRUCTURAL SCAFFOLD ONLY.
//
// This file defines the public API and data structures for an F4 engine that
// would drop into the existing rl_server / strategy framework. The hot inner
// pieces (symbolic preprocessing, Macaulay matrix construction, sparse
// reduction) are stubbed with `unimplemented!` and detailed TODOs.
//
// Why F4 vs Buchberger:
//   - Buchberger processes one S-polynomial at a time.
//   - F4 batches ALL pairs of a given lcm-degree into one matrix and reduces
//     them simultaneously via Gaussian elimination. Same Gröbner basis output,
//     ~10× faster on cyclic-n and dramatically faster on Semaev systems with
//     m ≥ 3 (where Buchberger gets stuck on basis-size blow-up).
//
// References for an implementor:
//   - Faugère (1999), "A new efficient algorithm for computing Gröbner bases (F4)".
//     The original. Sections 2-3 give Macaulay matrix construction; section 4
//     covers symbolic preprocessing.
//   - msolve (https://msolve.lip6.fr/) — open-source C implementation. Source
//     is the cleanest way to see modern engineering tradeoffs (sparse linbox,
//     finite-field SIMD).
//   - openf4 (https://github.com/nauotit/openf4) — smaller reference impl.

use crate::monomial::Monomial;
use crate::pair::CriticalPair;
use crate::poly::Poly;
use crate::spoly::spoly;

// Mirrors BuchbergerState's field layout so the RL env can swap engines with
// no env.rs changes.
pub struct F4State {
    pub basis: Vec<Poly>,
    pub sugars: Vec<u32>,
    pub pairs: Vec<CriticalPair>,
    pub nvars: usize,
    pub step_count: usize,
    pub arith_ops: usize,
}

pub struct F4StepResult {
    /// Number of new basis polynomials added by this matrix-reduction round.
    pub added: usize,
    pub matrix_rows: usize,
    pub matrix_cols: usize,
}

impl F4State {
    pub fn new(initial: Vec<Poly>) -> Self {
        let nvars = initial.first().map(|p| p.nvars).unwrap_or(0);
        let sugars: Vec<u32> = initial.iter().map(|p| p.total_degree()).collect();
        let mut state = F4State {
            basis: initial,
            sugars,
            pairs: Vec::new(),
            nvars,
            step_count: 0,
            arith_ops: 0,
        };
        // Identical to BuchbergerState::new initial pair gen.
        for i in 0..state.basis.len() {
            for j in (i + 1)..state.basis.len() {
                state.install_pair_raw(i, j);
            }
        }
        state
    }

    fn install_pair_raw(&mut self, i: usize, j: usize) {
        // Same coprime-LM short-circuit as Buchberger; full GM is later.
        let (lmi, lmj) = match (self.basis[i].lm(), self.basis[j].lm()) {
            (Some(a), Some(b)) => (a.clone(), b.clone()),
            _ => return,
        };
        if lmi.gcd_is_one(&lmj) { return; }
        let lcm = lmi.lcm(&lmj);
        let dfi = lcm.degree() - lmi.degree();
        let dfj = lcm.degree() - lmj.degree();
        let sugar = (self.sugars[i] + dfi).max(self.sugars[j] + dfj);
        self.pairs.push(CriticalPair { i, j, lcm, sugar });
    }

    pub fn is_done(&self) -> bool { self.pairs.is_empty() }

    /// Process ONE F4 round: pick all pairs at the minimum lcm degree, build
    /// the Macaulay matrix, reduce, install nonzero rows as new basis polys.
    ///
    /// The RL action analog is: which *degree bucket* to drain next. (Plain F4
    /// always picks the minimum; signature-based variants like F5 use the
    /// Möller-Mora ordering. A learned policy could pick across non-min degrees.)
    pub fn step(&mut self, degree_idx: usize) -> F4StepResult {
        self.step_count += 1;
        let selected = self.select_pairs_by_degree(degree_idx);
        let s_polys = self.spolynomials_from_pairs(&selected);
        let (rows, columns) = symbolic_preprocess(s_polys, &self.basis);
        let reduced = matrix_reduce(rows, &columns);
        let mut added = 0;
        for p in reduced {
            if !p.is_zero() && !already_in_basis_as_lm(&p, &self.basis) {
                let idx = self.basis.len();
                self.basis.push(p);
                self.sugars.push(self.basis[idx].total_degree());
                gebauer_moller_install(self, idx);
                added += 1;
            }
        }
        F4StepResult { added, matrix_rows: 0, matrix_cols: 0 }
    }

    /// Bucket the queued pairs by lcm degree, return the indices in
    /// bucket # `degree_idx` (0 = lowest degree). The default F4 strategy is
    /// always degree_idx = 0 (drain the lowest first).
    fn select_pairs_by_degree(&mut self, _degree_idx: usize) -> Vec<CriticalPair> {
        // TODO: group self.pairs by lcm.degree(), pop the chosen bucket.
        Vec::new()
    }

    fn spolynomials_from_pairs(&self, pairs: &[CriticalPair]) -> Vec<Poly> {
        pairs.iter()
            .map(|p| spoly(&self.basis[p.i], &self.basis[p.j]))
            .collect()
    }
}

// === Stubs that need real implementation ==============================

/// Symbolic preprocessing (Faugère 1999, §3.3).
///
/// Given a set of polynomials to reduce against `basis`, walk each polynomial's
/// non-leading monomials, and for every monomial m that has a basis poly g with
/// LM(g) dividing m, add the *shifted* g (m/LM(g) * g) as an extra row so the
/// matrix-stage reduction can cancel m. Returns:
///   - rows: S-polys ∪ all shifted reductors
///   - columns: the union of monomials appearing in any row, sorted descending.
///
/// TODO: implement. Use a worklist over monomials to discover new reductors.
pub fn symbolic_preprocess(_s_polys: Vec<Poly>, _basis: &[Poly])
    -> (Vec<Poly>, Vec<Monomial>)
{
    unimplemented!("symbolic_preprocess: see Faugère 1999 §3.3")
}

/// Macaulay matrix reduction. `rows` are polynomials viewed over the global
/// monomial set `columns` (sorted descending). Performs row-echelon reduction
/// over Fp, returns the reduced rows (still as Polys).
///
/// TODO: dense implementation first (cubic, easy), then sparse. For Fp this is
/// just Gaussian elimination with pivots ordered by leading column. Build a
/// matrix `M[B][cols] : Fp` from the rows, run REF, extract back to polys.
/// Track arith_ops for parity with BuchbergerState reporting.
pub fn matrix_reduce(_rows: Vec<Poly>, _columns: &[Monomial]) -> Vec<Poly> {
    unimplemented!("matrix_reduce: row-echelon reduction of Macaulay matrix")
}

/// Quick check: does `p`'s leading monomial coincide with an existing basis
/// LM? (Same LM → already in the ideal up to a unit; don't re-add.)
fn already_in_basis_as_lm(p: &Poly, basis: &[Poly]) -> bool {
    let lm = match p.lm() { Some(m) => m, None => return false };
    basis.iter().any(|g| g.lm() == Some(lm))
}

/// Reuse the GM pruning logic from buchberger.rs when a new poly is added.
/// In a real impl this would be either copied in or refactored to a shared
/// helper. Stubbed here to keep this file self-contained.
fn gebauer_moller_install(_state: &mut F4State, _new_idx: usize) {
    // TODO: identical to BuchbergerState::apply_gebauer_moller. Once F4
    // matures, lift GM into its own module shared between the two engines.
}

// === Notes on RL integration ==========================================
//
// Strategy hook in F4 is different from Buchberger: instead of "which pair",
// it's "which degree bucket" (and within a degree bucket, optionally "which
// symbolic-preprocess reductors to include"). The current pair-feature
// pipeline doesn't directly transfer. A future env.rs would expose:
//   - degree_buckets: Vec<(degree, n_pairs_in_bucket, avg_sugar, max_sugar)>
//   - and action ∈ {0..n_buckets}, defaulting to 0 = lowest degree.
// The PointerPolicy architecture (per-item embedding → logit) works as-is.

#[cfg(test)]
mod tests {
    // TODO: once symbolic_preprocess + matrix_reduce are implemented, port the
    // cyclic-n tests and verify F4 produces the same Gröbner basis (same LMs)
    // as BuchbergerState on cyclic-3..5.
}
