//! # Refutation-degree (last-fall / PC-degree) harness for binary Semaev PDP.
//!
//! Companion to [`crate::cryptanalysis::ffd_harness`].  Where that module
//! measures the **first** fall degree `d_ff` of the unrestricted binary
//! Semaev `S₃` system, this one measures the number that actually decides
//! security: the **refutation degree** `D*` — the smallest Macaulay
//! degree at which a Gröbner / Polynomial-Calculus computation derives the
//! unit `1` from a **non-decomposable** point-decomposition instance.
//!
//! ## Why a different instance
//!
//! `ffd_harness` lets `X₁, X₂` range over *all* of `F_{2ⁿ}` (the full `2n`
//! bit-variables).  That system is generically *satisfiable* — `S₃ = 0`
//! defines a curve with ≈ `2ⁿ` points — so `1 ∉ ⟨S₃, field eqs⟩` and there
//! is **no refutation** to measure.  The genuine index-calculus point-
//! decomposition problem (PDP) restricts the factor base to a proper
//! `F_2`-subspace `V ⊂ F_{2ⁿ}` of dimension `n' < n` and asks
//!
//! ```text
//!     ∃ X₁, X₂ ∈ V :  S₃(X₁, X₂, x(R)) = 0 ?
//! ```
//!
//! For `2n' < n` *most* targets `x(R)` are **non-decomposable**: the
//! restricted system is *unsatisfiable*, so `1 ∈ ⟨ … ⟩` and a Gröbner
//! computation must derive it.  By the Clegg–Edmonds–Impagliazzo
//! correspondence the degree at which `1` first appears in the Macaulay
//! row-space is exactly the Polynomial-Calculus refutation degree, which
//! equals (up to the usual `O(1)`) the Huang–Kosters–Yeo **last fall
//! degree** and the Gröbner solving degree.  See
//! `RESEARCH_FFD_PROOF_COMPLEXITY.md` §2.
//!
//! ## What this measures, operationally
//!
//! Per `(n, n')` and a non-decomposable target `x₃`:
//!
//! 1. Build the `n` quadratic `F_2`-equations via
//!    [`ffd_harness::weil_descend_s3`], then **restrict** to the subspace
//!    `V = span{z⁰, …, z^{n'-1}}` by forcing the high bits of `X₁, X₂` to
//!    zero ([`restrict_to_subspace`]) — leaving `2n'` bit-variables.
//! 2. At each degree `D = 2, 3, …`, build the Macaulay matrix
//!    ([`ffd_harness::build_macaulay_rows`]) and test whether the constant
//!    monomial `1` (column 0) lies in the row-space ([`e0_in_rowspace`]).
//! 3. Report:
//!    - **first fall degree** `d_ff` — smallest `D` with `rank < rows` and
//!      `rank < cols` (same operational definition as `ffd_harness`), and
//!    - **refutation degree** `D*` — smallest `D` with `1 ∈ row-space`.
//!
//! The gap `D* − d_ff` is the Huang–Kosters–Yeo invariant made
//! operational, and the proposal's **prediction #1** is that across a
//! field sweep `d_ff` stays ≈ constant while `D*` climbs.
//!
//! ## Calibration
//!
//! Decomposability of `x₃` over `V` is checked by brute force over
//! `V × V` ([`is_decomposable`]).  A decomposable target is *satisfiable*,
//! so it must show **no** refutation (`D* = None`); a non-decomposable
//! target must eventually refute.  The test
//! `refutation_iff_nondecomposable` asserts exactly this equivalence at
//! small `n`, which validates the whole pipeline.
//!
//! ## Limits
//!
//! - Cost is dominated by the Macaulay rank at degree `D` in `v = 2n'`
//!   variables.  Comfortable to `v ≤ 10` (`n ≤ 10` with `n' = ⌊n/2⌋`) and
//!   `D ≤ v`.  Beyond that, swap in a sparse F4 backend
//!   (`cryptanalysis::groebner_f4`).
//! - Single `(b, x₃)` per `(n, n')`; for paper-grade numbers average over
//!   many non-decomposable targets.
//!
//! ## References
//!
//! - M. Clegg, J. Edmonds, R. Impagliazzo, *Using the Groebner basis
//!   algorithm to find proofs of unsatisfiability*, STOC 1996.
//! - M.-D. Huang, M. Kosters, S. L. Yeo, *Last fall degree, HFE, and Weil
//!   descent attacks on ECDLP*, CRYPTO 2015.
//! - I. Semaev, *Summation polynomials …*, ePrint 2004/031.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::binary_semaev::binary_semaev_s3;
use crate::cryptanalysis::ffd_harness::{
    build_macaulay_rows, choose_irreducible, f2_rank, quad_monomial_index, random_nonzero_f2m,
    weil_descend_s3, F2BoolPoly,
};
use rand::rngs::StdRng;
use rand::SeedableRng;

// ── Public API ──────────────────────────────────────────────────────

/// One row of the refutation-degree table.
#[derive(Clone, Debug)]
pub struct PcRow {
    /// Field extension degree `n` (system lives in `F_{2ⁿ}`).
    pub n: u32,
    /// Factor-base subspace dimension `n'` (`V = span{z⁰..z^{n'-1}}`).
    pub n_sub: u32,
    /// Bit-variables after restriction (`2n'`).
    pub num_vars: u32,
    /// `F_2`-equations after Weil descent (`n`).
    pub num_eqs: u32,
    /// Whether the target `x₃` decomposes over `V` (brute-force checked).
    /// Non-decomposable ⇒ the system is unsatisfiable ⇒ a refutation
    /// exists.
    pub decomposable: bool,
    /// Per-degree Macaulay measurements.
    pub per_degree: Vec<PcMeasurement>,
    /// Smallest `D` with a non-trivial fall (`rank < rows`, `rank < cols`)
    /// — the operational first-fall degree, comparable to `ffd_harness`.
    pub first_fall: Option<u32>,
    /// Smallest `D` at which the constant `1` enters the row-space — the
    /// refutation degree `D*` (≈ last fall degree / solving degree).
    /// `None` if no refutation up to `d_max` (expected when `decomposable`).
    pub refutation_degree: Option<u32>,
}

#[derive(Clone, Debug)]
pub struct PcMeasurement {
    pub degree: u32,
    pub rows_constructed: u64,
    pub cols: u64,
    pub rank: u64,
    /// `true` iff the constant monomial `1` is in the Macaulay row-space
    /// at this degree (i.e. a degree-`D` refutation exists).
    pub refuted: bool,
}

/// Measure the refutation degree for a single `(n, n', b, x₃)` instance.
pub fn measure_refutation_degree(
    n: u32,
    n_sub: u32,
    irr: &IrreduciblePoly,
    b: &F2mElement,
    x3: &F2mElement,
    d_max: u32,
) -> PcRow {
    assert!(n_sub >= 1 && n_sub <= n, "need 1 ≤ n' ≤ n");
    let full_eqs = weil_descend_s3(n, irr, b, x3);
    let eqs = restrict_to_subspace(&full_eqs, n, n_sub);
    let num_vars = 2 * n_sub;
    let num_eqs = eqs.len() as u32;
    let decomposable = is_decomposable(n, n_sub, irr, b, x3);

    let (first_fall, refutation_degree, per_degree) = refutation_scan(&eqs, num_vars, d_max);

    PcRow {
        n,
        n_sub,
        num_vars,
        num_eqs,
        decomposable,
        per_degree,
        first_fall,
        refutation_degree,
    }
}

/// Scan Macaulay degrees `2..=d_max` of an arbitrary `F_2`-quadratic system
/// `eqs` in `num_vars` Boolean variables, returning
/// `(first_fall, refutation_degree, per_degree)`:
///
/// - `first_fall` = smallest `D` with `rank < rows` and `rank < cols`
///   (operational first-fall onset), and
/// - `refutation_degree` = smallest `D` at which the constant `1` enters
///   the row-space (the PC / last-fall refutation degree `D*`).
///
/// Factored out of [`measure_refutation_degree`] so callers that build the
/// descended system differently — e.g. EXP-E's arbitrary factor-base
/// subspace via [`crate::cryptanalysis::descent_lowgamma`] — can reuse the
/// exact same measurement.
pub fn refutation_scan(
    eqs: &[F2BoolPoly],
    num_vars: u32,
    d_max: u32,
) -> (Option<u32>, Option<u32>, Vec<PcMeasurement>) {
    let mut per_degree = Vec::new();
    let mut first_fall: Option<u32> = None;
    let mut refutation_degree: Option<u32> = None;

    for d in 2..=d_max {
        let (rows, cols, rows_constructed) = build_macaulay_rows(eqs, num_vars, d);
        // Rank of the Macaulay matrix (clone — f2_rank reduces in place).
        let mut for_rank = rows.clone();
        let rank = f2_rank(&mut for_rank, cols) as u64;
        let refuted = e0_in_rowspace(&rows, cols, rank as usize);

        let nontrivial_syzygy = rank < rows_constructed;
        let not_saturated = rank < cols as u64;
        if first_fall.is_none() && nontrivial_syzygy && not_saturated {
            first_fall = Some(d);
        }
        if refutation_degree.is_none() && refuted {
            refutation_degree = Some(d);
        }

        per_degree.push(PcMeasurement {
            degree: d,
            rows_constructed,
            cols: cols as u64,
            rank,
            refuted,
        });

        // Once refuted, deeper degrees only re-confirm it — stop early.
        if refutation_degree.is_some() {
            break;
        }
    }
    (first_fall, refutation_degree, per_degree)
}

/// Run the refutation-degree sweep over `n_range`.  For each `n` the
/// factor-base dimension is `n' = ⌊n/2⌋` (so `|V|² ≈ 2ⁿ` and roughly half
/// the targets are non-decomposable), and a **non-decomposable** target is
/// found by scanning random `x₃` (so every reported `D*` is a genuine
/// refutation).  `d_max` caps the Macaulay degree; `seed` makes the random
/// `(b, x₃)` choices reproducible.
pub fn run_pc_sweep(n_range: std::ops::RangeInclusive<u32>, d_max: u32, seed: u64) -> Vec<PcRow> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out = Vec::new();
    for n in n_range {
        let n_sub = (n / 2).max(1);
        let irr = choose_irreducible(n);
        let b = random_nonzero_f2m(&mut rng, n);
        // Find a non-decomposable target (bounded attempts; fall back to
        // the last candidate if the field is too small to have one).
        let mut x3 = random_nonzero_f2m(&mut rng, n);
        for _ in 0..256 {
            if !is_decomposable(n, n_sub, &irr, &b, &x3) {
                break;
            }
            x3 = random_nonzero_f2m(&mut rng, n);
        }
        out.push(measure_refutation_degree(n, n_sub, &irr, &b, &x3, d_max));
    }
    out
}

// ── Subspace restriction ────────────────────────────────────────────

/// Restrict a Weil-descent system over the full `2n` bit-variables to the
/// factor-base subspace `V = span{z⁰, …, z^{n'-1}}` by setting the high
/// bits of `X₁` and `X₂` (positions `n'..n`) to zero, then re-indexing the
/// surviving `2n'` variables as `[X₁ low bits ‖ X₂ low bits]`.
///
/// This is a literal `x_i ← 0` substitution: any monomial that touches a
/// dropped variable vanishes.
pub fn restrict_to_subspace(eqs: &[F2BoolPoly], n: u32, n_sub: u32) -> Vec<F2BoolPoly> {
    let old_vars = 2 * n;
    let new_vars = 2 * n_sub;
    eqs.iter()
        .map(|eq| {
            let mut out = F2BoolPoly::zero(new_vars);
            // Constant term.
            out.coeffs[0] ^= eq.coeffs[0];
            // Linear terms.
            for i in 0..old_vars {
                if eq.coeffs[1 + i as usize] {
                    if let Some(m) = map_var(i, n, n_sub) {
                        out.coeffs[1 + m as usize] ^= true;
                    }
                }
            }
            // Quadratic terms.
            for i in 0..old_vars {
                for j in (i + 1)..old_vars {
                    let idx = quad_monomial_index(i, j, old_vars);
                    if idx < eq.coeffs.len() && eq.coeffs[idx] {
                        if let (Some(mi), Some(mj)) =
                            (map_var(i, n, n_sub), map_var(j, n, n_sub))
                        {
                            // mi != mj because the var map is injective.
                            let (a, c) = if mi < mj { (mi, mj) } else { (mj, mi) };
                            let nidx = quad_monomial_index(a, c, new_vars);
                            out.coeffs[nidx] ^= true;
                        }
                    }
                }
            }
            out
        })
        .collect()
}

/// Map an old bit-variable index (in the `2n`-variable layout
/// `[X₁ bits 0..n ‖ X₂ bits 0..n]`) to its index in the restricted
/// `2n'`-variable layout, or `None` if it is a dropped high bit.
fn map_var(old: u32, n: u32, n_sub: u32) -> Option<u32> {
    if old < n_sub {
        Some(old) // X₁ low bit
    } else if old < n {
        None // X₁ high bit → forced 0
    } else if old < n + n_sub {
        Some(n_sub + (old - n)) // X₂ low bit
    } else {
        None // X₂ high bit → forced 0
    }
}

// ── Decomposability oracle (brute force over V × V) ─────────────────

/// Does `x₃` admit a decomposition `S₃(X₁, X₂, x₃) = 0` with both
/// `X₁, X₂ ∈ V = span{z⁰, …, z^{n'-1}}`?  Brute force over `|V|² = 2^{2n'}`
/// pairs — cheap for the `n' ≤ 6` regime this harness targets.
pub fn is_decomposable(
    n: u32,
    n_sub: u32,
    irr: &IrreduciblePoly,
    b: &F2mElement,
    x3: &F2mElement,
) -> bool {
    let span = 1u32 << n_sub;
    for mask1 in 0..span {
        let x1 = subspace_element(mask1, n);
        for mask2 in 0..span {
            let x2 = subspace_element(mask2, n);
            if binary_semaev_s3(&x1, &x2, x3, b, irr).is_zero() {
                return true;
            }
        }
    }
    false
}

/// Build the `F_{2ⁿ}` element whose low `n'` bits are given by `mask`
/// (bit `k` of `mask` sets coordinate `z^k`).
fn subspace_element(mask: u32, n: u32) -> F2mElement {
    let positions: Vec<u32> = (0..32).filter(|k| (mask >> k) & 1 == 1).collect();
    F2mElement::from_bit_positions(&positions, n)
}

// ── Refutation test: is the unit `1` in the Macaulay row-space? ─────

/// Test whether the constant monomial `1` (column 0) lies in the row-space
/// of the Macaulay matrix.  Equivalent to: appending the unit row `e₀`
/// does **not** increase the rank.  `base_rank` is the already-computed
/// rank of `rows`, so this costs one extra rank reduction.
pub fn e0_in_rowspace(rows: &[Vec<u64>], cols: usize, base_rank: usize) -> bool {
    if rows.is_empty() || cols == 0 {
        return false;
    }
    let words = (cols + 63) / 64;
    let mut augmented = rows.to_vec();
    let mut e0 = vec![0u64; words];
    e0[0] = 1; // column 0 = the constant monomial `1`
    augmented.push(e0);
    let aug_rank = f2_rank(&mut augmented, cols);
    aug_rank == base_rank
}

// ── Pretty printing ─────────────────────────────────────────────────

pub fn print_pc_sweep(rows: &[PcRow]) {
    println!();
    println!(
        "{:>3} {:>3} {:>5} {:>4} {:>6} | {:>9} {:>9}",
        "n", "n'", "vars", "eqs", "decmp", "first-fall", "D* (refute)"
    );
    println!("{}", "─".repeat(58));
    for r in rows {
        let ff = r
            .first_fall
            .map(|d| d.to_string())
            .unwrap_or_else(|| "—".into());
        let ds = r
            .refutation_degree
            .map(|d| d.to_string())
            .unwrap_or_else(|| "none".into());
        println!(
            "{:>3} {:>3} {:>5} {:>4} {:>6} | {:>9} {:>9}",
            r.n,
            r.n_sub,
            r.num_vars,
            r.num_eqs,
            if r.decomposable { "yes" } else { "no" },
            ff,
            ds
        );
    }
    println!();
    println!(
        "Reading: `first-fall` is the syzygy onset measured by ffd_harness;\n\
         `D*` is the Polynomial-Calculus refutation degree (≈ last fall\n\
         degree). Prediction #1: D* climbs with n while first-fall stays\n\
         roughly flat. All `D*` rows are non-decomposable (genuinely\n\
         unsatisfiable) instances."
    );
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn irr(n: u32) -> IrreduciblePoly {
        choose_irreducible(n)
    }

    /// Restriction drops every monomial that touches a high bit and
    /// re-indexes into `2n'` variables.
    #[test]
    fn restriction_drops_high_bit_monomials() {
        let n = 4;
        let n_sub = 2;
        let irr = irr(n);
        let b = F2mElement::from_bit_positions(&[0, 1], n);
        let x3 = F2mElement::from_bit_positions(&[1, 2], n);
        let full = weil_descend_s3(n, &irr, &b, &x3);
        let restricted = restrict_to_subspace(&full, n, n_sub);
        // Same number of equations, smaller variable space.
        assert_eq!(restricted.len(), full.len());
        let new_len = F2BoolPoly::zero(2 * n_sub).coeffs.len();
        for eq in &restricted {
            assert_eq!(eq.coeffs.len(), new_len);
        }
    }

    /// `map_var` is injective on the kept indices and `None` on dropped
    /// high bits.
    #[test]
    fn map_var_layout() {
        let (n, n_sub) = (5, 2);
        // X₁ low bits 0,1 → 0,1 ; X₂ low bits 5,6 → 2,3.
        assert_eq!(map_var(0, n, n_sub), Some(0));
        assert_eq!(map_var(1, n, n_sub), Some(1));
        assert_eq!(map_var(2, n, n_sub), None); // X₁ high bit
        assert_eq!(map_var(4, n, n_sub), None);
        assert_eq!(map_var(5, n, n_sub), Some(2)); // X₂ bit 0
        assert_eq!(map_var(6, n, n_sub), Some(3));
        assert_eq!(map_var(7, n, n_sub), None); // X₂ high bit
    }

    /// **The calibration theorem**: over a small field, a target refutes
    /// (`1` enters the row-space) **iff** it is non-decomposable over `V`.
    /// This validates that `D*` measures the right thing.
    #[test]
    fn refutation_iff_nondecomposable() {
        let n = 4;
        let n_sub = 2; // |V| = 4, |V|² = 16 pairs
        let irr = irr(n);
        let b = F2mElement::from_bit_positions(&[0], n);
        let d_max = 2 * n_sub + 1; // enough degree to saturate at n'=2

        let mut checked = 0;
        for t in 1u32..(1u32 << n) {
            let x3 = subspace_element_full(t, n);
            if x3.is_zero() {
                continue;
            }
            let row = measure_refutation_degree(n, n_sub, &irr, &b, &x3, d_max);
            // The equivalence: decomposable ⇔ NOT refuted.
            assert_eq!(
                row.decomposable,
                row.refutation_degree.is_none(),
                "x3=#{t}: decomposable={} but refutation={:?}",
                row.decomposable,
                row.refutation_degree
            );
            checked += 1;
        }
        assert!(checked > 0, "no targets were checked");
    }

    /// Build an `F_{2ⁿ}` element from the low bits of `t` (all `n` bits,
    /// not just the subspace) — used to range over every target.
    fn subspace_element_full(t: u32, n: u32) -> F2mElement {
        let positions: Vec<u32> = (0..n).filter(|k| (t >> k) & 1 == 1).collect();
        F2mElement::from_bit_positions(&positions, n)
    }

    /// A non-decomposable target produces a finite refutation degree
    /// `D* ≥ 2`, and the decomposability oracle agrees.
    #[test]
    fn nondecomposable_target_refutes() {
        let n = 4;
        let n_sub = 2;
        let irr = irr(n);
        let b = F2mElement::from_bit_positions(&[0], n);
        // Search for a non-decomposable target.
        let mut found = None;
        for t in 1u32..(1u32 << n) {
            let x3 = subspace_element_full(t, n);
            if x3.is_zero() {
                continue;
            }
            if !is_decomposable(n, n_sub, &irr, &b, &x3) {
                found = Some(x3);
                break;
            }
        }
        let x3 = found.expect("expected a non-decomposable target at n=4, n'=2");
        let row = measure_refutation_degree(n, n_sub, &irr, &b, &x3, 2 * n_sub + 1);
        assert!(!row.decomposable);
        let dstar = row.refutation_degree.expect("non-decomposable must refute");
        assert!(dstar >= 2);
    }

    /// `e0_in_rowspace` is correct on tiny hand-built matrices.
    #[test]
    fn e0_membership_basic() {
        // cols = 3. Rows {110, 011}: span does NOT contain e0 = 100
        // (every combination has even-or-paired structure: 110, 011,
        // 110^011 = 101 — none is 100).
        let rows = vec![vec![0b110u64], vec![0b011u64]];
        let base = f2_rank(&mut rows.clone(), 3);
        assert!(!e0_in_rowspace(&rows, 3, base));
        // Add row 010: now 110 ^ 010 = 100 = e0 is in the span.
        let rows2 = vec![vec![0b110u64], vec![0b011u64], vec![0b010u64]];
        let base2 = f2_rank(&mut rows2.clone(), 3);
        assert!(e0_in_rowspace(&rows2, 3, base2));
    }

    /// End-to-end sweep runs and every reported `D*` is from a
    /// non-decomposable (unsatisfiable) instance.
    #[test]
    fn pc_sweep_runs_and_is_nondecomposable() {
        let rows = run_pc_sweep(4..=6, 7, 0x_DEC0_DE_u64);
        assert_eq!(rows.len(), 3);
        for r in &rows {
            assert_eq!(r.num_eqs, r.n);
            assert_eq!(r.num_vars, 2 * r.n_sub);
            if let Some(_d) = r.refutation_degree {
                // A genuine refutation requires an unsatisfiable instance.
                assert!(!r.decomposable, "refuted a decomposable target?!");
            }
        }
    }
}
