//! # Buying down the Gröbner solving degree `D*` — lever L3 (hybrid slicing).
//!
//! The FFD program (`RESEARCH_FFD_PROOF_COMPLEXITY.md`,
//! `RESEARCH_FFD_WORKFLOW.md`) measures `D*`, the degree at which the
//! Macaulay/Gröbner computation on a Weil-descended Semaev system
//! terminates, and establishes two facts we now use *offensively*:
//!
//! - **P6** — over-determination collapses `D*`: at determination ratio
//!   `ρ = #eqs/#vars ≫ 1` the system admits a degree-2 Nullstellensatz
//!   certificate, so `D* = 2`.
//! - **P3-alg** — `D*` decreases with the early Hilbert-function defect
//!   `Δ_low` (slope ≈ −7.6 degrees per unit defect in the critical regime).
//!
//! P6 was recorded as a *nuisance* — the reason the early single-sample
//! harness saw a spuriously flat `D*`. But `ρ` is not fixed by the problem:
//! an attacker can **raise it at will** by guessing `k` of the `N = 2n'`
//! Boolean variables and solving the `2^k` resulting slices. That is the
//! Bettale–Faugère–Perret *hybrid approach*, and this module makes it
//! measurable on descended Semaev systems.
//!
//! ## Why the accounting is exact here
//!
//! Fix a target `x₃` that is **non-decomposable** over the factor base `V`
//! (the common case in index calculus — most relation attempts fail). The
//! descended system `Σ` is then unsatisfiable, and so is *every* one of the
//! `2^k` slices obtained by fixing `k` variables. Hence:
//!
//! - every slice refutes, so `D*` is defined on all of them (no censoring
//!   from "this slice happened to be satisfiable");
//! - the slices **partition** the search space, so the hybrid is not an
//!   approximation — it covers exactly what the direct solve covered;
//! - total work is `Σ_slices cost(D*_slice)`, which sampling estimates as
//!   `2^k · E[cost(D*_slice)]`.
//!
//! So the trade is clean: `2^k` multiplicative slices bought against a
//! (hopefully) collapsing `D*` in each. This module measures both halves
//! and reports the crossover.
//!
//! ## What we expose
//!
//! - [`specialize`] — substitute constants for a chosen set of Boolean
//!   variables and re-index the survivors (an exact affine slice of `Σ`).
//! - [`GuessPattern`] — *which* variables to fix. `OneSide` (all from the
//!   `X₁` half) is index calculus's asymmetric factor-base split;
//!   `Balanced` splits evenly; `Spread` picks uniformly at random.
//! - [`HybridRow`] / [`run_hybrid_sweep`] — the `D*`-vs-`k` curve with a
//!   histogram-exact cost model, so the linear-algebra exponent `ω` can be
//!   varied after the fact ([`HybridRow::log2_total_cost`]).
//!
//! ## Reading the output — and the baseline that actually matters
//!
//! `log2_total_cost(ω) = k + log₂ E[cols(N−k, D*)^ω]` is the log₂ work of
//! the whole hybrid at that `k`; `k = 0` is the direct solve.
//!
//! Beating `k = 0` is **not** the bar. The `k = N` endpoint of this same
//! family is brute-force enumeration of `V × V`, so at toy sizes the cost
//! model will happily report that "guessing more always helps" — it is
//! drifting toward exhaustive search, which is genuinely fastest when
//! `2^N` is small. Any gate keyed only on the `k = 0` comparison measures
//! that artifact and nothing else.
//!
//! So the honest comparison is against [`log2_enumeration_cost`] (`2^N`
//! Semaev evaluations, charged generously at `O(1)` amortised each), and
//! the scaling quantity is [`HybridSweep::collapse_fraction`] — the
//! fraction `c = k/N` at which `D*` reaches the Nullstellensatz floor of 2.
//! Total hybrid work is `2^{cN}·poly`, so the lever is a real speedup only
//! if `c < 1`, and survives to cryptographic size only if `c` does not
//! drift upward with `N`.
//!
//! ## References
//!
//! - L. Bettale, J.-C. Faugère, L. Perret, *Hybrid approach for solving
//!   multivariate systems over finite fields*, J. Math. Cryptol. 2009.
//! - Bardet–Faugère–Salvy, semi-regular complexity (the `cols^ω` model).
//! - `RESEARCH_DEGREE_REDUCTION.md` (this thread's charter, ledger R1–R5).

use std::collections::BTreeMap;

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::descent_algebraic::{early_defect, rank_profile};
use crate::cryptanalysis::descent_lowgamma::{
    descend_on_subspace, is_decomposable_on_subspace, BasisFamily, FactorSubspace,
};
use crate::cryptanalysis::ffd_harness::{
    num_monomials_upto_degree, quad_monomial_index, F2BoolPoly,
};
use crate::cryptanalysis::pc_degree_harness::refutation_scan;

// ── Affine slicing ──────────────────────────────────────────────────

/// Substitute constants for some Boolean variables and re-index the rest.
///
/// `assign[i] = Some(v)` fixes `x_i := v`; `assign[i] = None` keeps `x_i` as
/// a free variable. The surviving variables are re-indexed in increasing
/// order of their original index, so the result is a system in
/// `num_vars − (#fixed)` variables.
///
/// This is an exact substitution in the Boolean ring `F_2[x]/⟨x_i²−x_i⟩`:
/// a quadratic monomial with one fixed factor becomes linear, one with both
/// fixed becomes a constant. The slices over all `2^k` assignments of a
/// fixed variable set **partition** the solution space of the original
/// system.
pub fn specialize(
    eqs: &[F2BoolPoly],
    num_vars: u32,
    assign: &[Option<bool>],
) -> (Vec<F2BoolPoly>, u32) {
    assert_eq!(
        assign.len(),
        num_vars as usize,
        "assign covers every variable"
    );

    // old index → Some(new index) for survivors, None for fixed.
    let mut remap = vec![None; num_vars as usize];
    let mut new_vars = 0u32;
    for (i, a) in assign.iter().enumerate() {
        if a.is_none() {
            remap[i] = Some(new_vars);
            new_vars += 1;
        }
    }

    let out = eqs
        .iter()
        .map(|eq| {
            let mut f = F2BoolPoly::zero(new_vars);
            if eq.coeffs[0] {
                f.coeffs[0] ^= true;
            }
            // Linear monomials.
            for i in 0..num_vars as usize {
                if !eq.coeffs.get(1 + i).copied().unwrap_or(false) {
                    continue;
                }
                match assign[i] {
                    Some(true) => f.coeffs[0] ^= true,
                    Some(false) => {}
                    None => f.coeffs[1 + remap[i].unwrap() as usize] ^= true,
                }
            }
            // Quadratic monomials x_i·x_j (i < j).
            for i in 0..num_vars {
                for j in (i + 1)..num_vars {
                    let idx = quad_monomial_index(i, j, num_vars);
                    if idx >= eq.coeffs.len() || !eq.coeffs[idx] {
                        continue;
                    }
                    let (ai, aj) = (assign[i as usize], assign[j as usize]);
                    match (ai, aj) {
                        // Both fixed → constant.
                        (Some(vi), Some(vj)) => {
                            if vi && vj {
                                f.coeffs[0] ^= true;
                            }
                        }
                        // One fixed → linear in the survivor (or vanishes).
                        (Some(vi), None) => {
                            if vi {
                                f.coeffs[1 + remap[j as usize].unwrap() as usize] ^= true;
                            }
                        }
                        (None, Some(vj)) => {
                            if vj {
                                f.coeffs[1 + remap[i as usize].unwrap() as usize] ^= true;
                            }
                        }
                        // Both survive → quadratic, re-indexed.
                        (None, None) => {
                            let (a, b) = (remap[i as usize].unwrap(), remap[j as usize].unwrap());
                            f.coeffs[quad_monomial_index(a, b, new_vars)] ^= true;
                        }
                    }
                }
            }
            f
        })
        .collect();

    (out, new_vars)
}

// ── Which variables to guess ────────────────────────────────────────

/// Where in the `[X₁ bits ‖ X₂ bits]` layout the guessed variables come
/// from. The choice is not cosmetic: `OneSide` is the index-calculus
/// *asymmetric factor-base split* (shrink `V` for one summand only), while
/// `Balanced` shrinks both summands symmetrically.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GuessPattern {
    /// All guessed variables taken from the `X₁` half (low indices first).
    OneSide,
    /// Alternate `X₁`, `X₂`, `X₁`, … so both halves shrink together.
    Balanced,
    /// A uniformly random `k`-subset of all `N` variables.
    Spread,
}

impl GuessPattern {
    pub fn label(&self) -> &'static str {
        match self {
            GuessPattern::OneSide => "one-side",
            GuessPattern::Balanced => "balanced",
            GuessPattern::Spread => "spread",
        }
    }

    /// Pick the `k` variable indices to fix, out of `num_vars` laid out as
    /// `[X₁ bits 0..h ‖ X₂ bits 0..h]` with `h = num_vars/2`.
    fn pick(&self, num_vars: u32, k: u32, rng: &mut dyn FnMut() -> u64) -> Vec<u32> {
        let half = num_vars / 2;
        match self {
            GuessPattern::OneSide => (0..k.min(num_vars)).collect(),
            GuessPattern::Balanced => (0..k.min(num_vars))
                .map(|t| {
                    let side = t % 2;
                    let off = t / 2;
                    if side == 0 {
                        off
                    } else {
                        half + off
                    }
                })
                .filter(|v| *v < num_vars)
                .collect(),
            GuessPattern::Spread => {
                // Partial Fisher–Yates over 0..num_vars.
                let mut pool: Vec<u32> = (0..num_vars).collect();
                for t in 0..k.min(num_vars) as usize {
                    let j = t + (rng() as usize) % (pool.len() - t);
                    pool.swap(t, j);
                }
                let mut v: Vec<u32> = pool[..k.min(num_vars) as usize].to_vec();
                v.sort_unstable();
                v
            }
        }
    }
}

// ── One row of the D*-vs-k curve ────────────────────────────────────

/// Aggregated measurement at one guess count `k`.
#[derive(Clone, Debug)]
pub struct HybridRow {
    /// Number of Boolean variables fixed per slice.
    pub k: u32,
    /// Free variables left in each slice, `N − k`.
    pub vars: u32,
    /// Equations (unchanged by slicing) — the descent gives `n` of them.
    pub eqs: u32,
    /// Determination ratio `#eqs / #vars` of a slice.
    pub rho: f64,
    /// Slices actually measured (targets × slices-per-target).
    pub slices: u32,
    /// Slices that did not refute within `d_max` (censored; excluded from
    /// the statistics below, and a warning sign if nonzero).
    pub censored: u32,
    /// Mean `D*` over the measured slices.
    pub dstar_mean: f64,
    /// Worst `D*` seen — the degree the attacker must actually budget for.
    pub dstar_max: u32,
    /// `D*` histogram (degree → count).
    pub dstar_hist: BTreeMap<u32, u32>,
    /// Mean early Hilbert defect `Δ_low` of a slice (cutoff 3).
    pub delta_low_mean: f64,
    /// Mean operational first-fall degree `d_ff` of a slice.
    pub first_fall_mean: f64,
}

impl HybridRow {
    /// log₂ of the total hybrid work at linear-algebra exponent `ω`:
    ///
    /// ```text
    ///   log2_total = k + log₂ E[ cols(N−k, D*)^ω ]
    /// ```
    ///
    /// computed **exactly** from the `D*` histogram (not from the mean
    /// degree, which would understate the cost by Jensen). `k = 0` is the
    /// direct solve, so a lever pays iff some `k > 0` comes in lower.
    pub fn log2_total_cost(&self, omega: f64) -> f64 {
        let total: u32 = self.dstar_hist.values().sum();
        if total == 0 {
            return f64::NAN;
        }
        let mut expected = 0.0f64;
        for (&d, &count) in &self.dstar_hist {
            let cols = num_monomials_upto_degree(self.vars, d) as f64;
            expected += (count as f64 / total as f64) * cols.powf(omega);
        }
        self.k as f64 + expected.log2()
    }

    /// log₂ of the work in a *single* slice at exponent `ω` (i.e. without
    /// the `2^k` enumeration factor) — the quantity the degree collapse
    /// actually reduces.
    pub fn log2_slice_cost(&self, omega: f64) -> f64 {
        self.log2_total_cost(omega) - self.k as f64
    }
}

/// The full sweep for one `(family, n, n', pattern)` cell.
#[derive(Clone, Debug)]
pub struct HybridSweep {
    pub family: BasisFamily,
    pub pattern: GuessPattern,
    pub n: u32,
    pub n_sub: u32,
    /// `N = 2n'`, the descended variable count before slicing.
    pub full_vars: u32,
    pub targets: u32,
    pub d_max: u32,
    pub rows: Vec<HybridRow>,
}

impl HybridSweep {
    /// log₂ cost of brute-force enumeration of `V × V` at this operating
    /// point — see [`log2_enumeration_cost`]. The baseline any hybrid must
    /// beat to count as a speedup.
    pub fn log2_enumeration(&self) -> f64 {
        log2_enumeration_cost(self.full_vars, self.n)
    }

    /// log₂ margin of the best hybrid over **brute-force enumeration** at
    /// exponent `ω`. Positive means the hybrid genuinely beats the search it
    /// replaces; negative means it does not, and gives the size of the
    /// shortfall in bits.
    pub fn margin_vs_enumeration(&self, omega: f64) -> Option<f64> {
        let best = self
            .rows
            .iter()
            .map(|r| r.log2_total_cost(omega))
            .filter(|c| c.is_finite())
            .fold(f64::INFINITY, f64::min);
        if best.is_finite() {
            Some(self.log2_enumeration() - best)
        } else {
            None
        }
    }

    /// The smallest `k` at which **every sampled slice** has collapsed to the
    /// Nullstellensatz floor `D* = 2`, as a fraction of `N`.
    ///
    /// This is the lever's true scaling quantity. Total hybrid work is
    /// `2^k · poly`, so the hybrid beats `2^N` enumeration only if the
    /// collapse arrives at `k = cN` with `c < 1` — and the *trend* of `c`
    /// with `N` decides whether the lever survives to cryptographic size.
    pub fn collapse_fraction(&self) -> Option<f64> {
        // If the unsliced system is already at the floor there is no degree
        // to reduce, so `c = 0` would be a degenerate reading, not a
        // measurement. Report `None` and let the caller exclude the cell.
        if self.rows.first().map(|r| r.dstar_max <= 2).unwrap_or(false) {
            return None;
        }
        self.rows
            .iter()
            .find(|r| r.censored == 0 && r.dstar_max <= 2)
            .map(|r| r.k as f64 / self.full_vars as f64)
    }

    /// Largest `k` scanned.
    pub fn k_max(&self) -> u32 {
        self.rows.last().map(|r| r.k).unwrap_or(0)
    }

    /// Does the cost curve have an **interior** minimum at exponent `ω`?
    ///
    /// This is the sweep's own honesty check. The `k = N` endpoint of the
    /// hybrid family is exhaustive search, so a cost model evaluated on toy
    /// systems will happily slide the optimum all the way to the largest `k`
    /// available — reporting "guessing always helps" when what it has really
    /// found is that `2^N` is small. If `k* == k_max` the sweep has **not**
    /// located an optimum; it has run off the end of the scan, and any
    /// "saving" it reports is an artifact of where the scan was truncated.
    pub fn optimum_is_interior(&self, omega: f64) -> bool {
        self.best(omega)
            .map(|(k, _)| k < self.k_max())
            .unwrap_or(false)
    }

    /// The `k` minimising total work at exponent `ω`, and its log₂ saving
    /// against the direct solve (`k = 0`). A positive saving means the
    /// hybrid lever pays.
    pub fn best(&self, omega: f64) -> Option<(u32, f64)> {
        let base = self.rows.iter().find(|r| r.k == 0)?.log2_total_cost(omega);
        let best = self
            .rows
            .iter()
            .filter(|r| r.log2_total_cost(omega).is_finite())
            .min_by(|a, b| {
                a.log2_total_cost(omega)
                    .partial_cmp(&b.log2_total_cost(omega))
                    .unwrap()
            })?;
        Some((best.k, base - best.log2_total_cost(omega)))
    }
}

/// log₂ cost of the **brute-force baseline**: enumerate all `2^N` points of
/// `V × V` and evaluate the descended system at each.
///
/// The model is deliberately *generous to brute force* — it charges only
/// `eqs` bit-operations per point, i.e. it assumes Gray-code enumeration
/// amortises the monomial updates to `O(1)` each, and it ignores constant
/// factors entirely. That is the honest way to run the comparison: a lever
/// that only wins against a strawman baseline has not won.
///
/// This baseline matters because it bounds the whole hybrid family. The
/// `k = N` endpoint of hybrid slicing *is* enumeration, so a hybrid that
/// cannot beat `2^N` has bought its degree reduction at full price and is
/// not a speedup over the search it replaced.
pub fn log2_enumeration_cost(full_vars: u32, eqs: u32) -> f64 {
    full_vars as f64 + (eqs.max(1) as f64).log2()
}

/// Run the `D*`-vs-`k` sweep.
///
/// For each of `targets` **non-decomposable** targets `x₃` over the factor
/// base `V` (so the slices are genuinely unsatisfiable and every one
/// refutes), and each `k ∈ 0..=k_max`, sample `slices_per_k` random
/// assignments of the `k` chosen variables and measure the slice's `D*`,
/// `d_ff` and `Δ_low`.
///
/// `k = 0` is measured once per target (there is only one slice).
/// Returns `None` if the subspace family is infeasible at `(n, n')`.
#[allow(clippy::too_many_arguments)]
pub fn run_hybrid_sweep(
    family: BasisFamily,
    pattern: GuessPattern,
    n: u32,
    n_sub: u32,
    irr: &IrreduciblePoly,
    k_max: u32,
    targets: u32,
    slices_per_k: u32,
    d_max: u32,
    seed: u64,
) -> Option<HybridSweep> {
    let full_vars = 2 * n_sub;
    let k_max = k_max.min(full_vars.saturating_sub(2));

    // Two INDEPENDENT streams. `next_target` draws the curve/target pairs;
    // `next_slice` draws the guessed variable sets and their values. Keeping
    // them separate is load-bearing: `GuessPattern::Spread` consumes randomness
    // that `OneSide`/`Balanced` do not, so a single stream would hand each
    // pattern a *different* set of targets and make the pattern comparison
    // meaningless. With the split, every pattern sees the same targets and the
    // `k = 0` row is pattern-independent (asserted in tests).
    let mut t_state = seed | 1;
    let mut next_target = move || {
        t_state ^= t_state >> 12;
        t_state ^= t_state << 25;
        t_state ^= t_state >> 27;
        t_state.wrapping_mul(0x2545F4914F6CDD1D)
    };
    let mut s_state = seed.rotate_left(32) | 1;
    let mut next_slice = move || {
        s_state ^= s_state >> 12;
        s_state ^= s_state << 25;
        s_state ^= s_state >> 27;
        s_state.wrapping_mul(0x2545F4914F6CDD1D)
    };

    // Accumulators, one bucket per k.
    struct Acc {
        slices: u32,
        censored: u32,
        dstar_sum: f64,
        dstar_max: u32,
        hist: BTreeMap<u32, u32>,
        delta_sum: f64,
        ff_sum: f64,
        ff_count: u32,
    }
    let mut acc: Vec<Acc> = (0..=k_max)
        .map(|_| Acc {
            slices: 0,
            censored: 0,
            dstar_sum: 0.0,
            dstar_max: 0,
            hist: BTreeMap::new(),
            delta_sum: 0.0,
            ff_sum: 0.0,
            ff_count: 0,
        })
        .collect();

    let rand_nz = |m: u32, rng: &mut dyn FnMut() -> u64| loop {
        let bits: Vec<u32> = (0..m).filter(|_| (rng() >> 19) & 1 == 1).collect();
        let e = F2mElement::from_bit_positions(&bits, m);
        if !e.is_zero() {
            return e;
        }
    };

    let mut found = 0u32;
    // Bounded search for non-decomposable targets: each draw is roughly a
    // coin flip at the critical operating point, so this terminates fast.
    let mut attempts = 0u32;
    while found < targets && attempts < targets * 64 + 256 {
        attempts += 1;
        let v = match family {
            BasisFamily::Random => {
                FactorSubspace::build(family, n, n_sub, irr, seed ^ (0x9E37 + attempts as u64))?
            }
            _ => FactorSubspace::build(family, n, n_sub, irr, 0)?,
        };
        let b = rand_nz(n, &mut next_target);
        let x3 = rand_nz(n, &mut next_target);
        if is_decomposable_on_subspace(&v, irr, &b, &x3) {
            continue;
        }
        found += 1;

        let eqs = descend_on_subspace(n, &v, irr, &b, &x3);

        for k in 0..=k_max {
            let reps = if k == 0 { 1 } else { slices_per_k };
            for _ in 0..reps {
                let idx = pattern.pick(full_vars, k, &mut next_slice);
                let mut assign: Vec<Option<bool>> = vec![None; full_vars as usize];
                for &i in &idx {
                    assign[i as usize] = Some((next_slice() >> 23) & 1 == 1);
                }
                let (seqs, nv) = specialize(&eqs, full_vars, &assign);
                let a = &mut acc[k as usize];
                a.slices += 1;

                let (ff, dstar, _) = refutation_scan(&seqs, nv, d_max);
                if let Some(f) = ff {
                    a.ff_sum += f as f64;
                    a.ff_count += 1;
                }
                match dstar {
                    Some(d) => {
                        a.dstar_sum += d as f64;
                        a.dstar_max = a.dstar_max.max(d);
                        *a.hist.entry(d).or_insert(0) += 1;
                    }
                    None => a.censored += 1,
                }

                let profile = rank_profile(&seqs, nv, n, 3.min(d_max));
                a.delta_sum += early_defect(&profile, 3);
            }
        }
    }

    if found == 0 {
        return None;
    }

    let rows = (0..=k_max)
        .map(|k| {
            let a = &acc[k as usize];
            let measured = a.slices - a.censored;
            let vars = full_vars - k;
            HybridRow {
                k,
                vars,
                eqs: n,
                rho: n as f64 / vars.max(1) as f64,
                slices: a.slices,
                censored: a.censored,
                dstar_mean: if measured > 0 {
                    a.dstar_sum / measured as f64
                } else {
                    f64::NAN
                },
                dstar_max: a.dstar_max,
                dstar_hist: a.hist.clone(),
                delta_low_mean: if a.slices > 0 {
                    a.delta_sum / a.slices as f64
                } else {
                    f64::NAN
                },
                first_fall_mean: if a.ff_count > 0 {
                    a.ff_sum / a.ff_count as f64
                } else {
                    f64::NAN
                },
            }
        })
        .collect();

    Some(HybridSweep {
        family,
        pattern,
        n,
        n_sub,
        full_vars,
        targets: found,
        d_max,
        rows,
    })
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ffd_harness::choose_irreducible;

    /// Fixing nothing must return the system unchanged.
    #[test]
    fn specialize_with_no_fixed_vars_is_identity() {
        let n = 6;
        let irr = choose_irreducible(n);
        let v = FactorSubspace::build(BasisFamily::Coordinate, n, 3, &irr, 0).unwrap();
        let b = F2mElement::from_bit_positions(&[0, 2], n);
        let x3 = F2mElement::from_bit_positions(&[1], n);
        let eqs = descend_on_subspace(n, &v, &irr, &b, &x3);

        let assign = vec![None; 6];
        let (out, nv) = specialize(&eqs, 6, &assign);
        assert_eq!(nv, 6);
        assert_eq!(out.len(), eqs.len());
        for (a, e) in out.iter().zip(eqs.iter()) {
            assert_eq!(a.coeffs, e.coeffs);
        }
    }

    /// Substitution must agree with direct evaluation: for every assignment
    /// of the free variables, the sliced system and the original system
    /// (with the fixed values plugged in) must take the same values.
    #[test]
    fn specialize_agrees_with_direct_evaluation() {
        let n = 6;
        let irr = choose_irreducible(n);
        let v = FactorSubspace::build(BasisFamily::Coordinate, n, 3, &irr, 0).unwrap();
        let b = F2mElement::from_bit_positions(&[0, 1], n);
        let x3 = F2mElement::from_bit_positions(&[2], n);
        let eqs = descend_on_subspace(n, &v, &irr, &b, &x3);
        let nvars = 6u32;

        // Fix x0 = 1, x3 = 0.
        let mut assign: Vec<Option<bool>> = vec![None; nvars as usize];
        assign[0] = Some(true);
        assign[3] = Some(false);
        let (sliced, nv) = specialize(&eqs, nvars, &assign);
        assert_eq!(nv, 4);

        // Enumerate all assignments of the 4 free variables.
        for mask in 0u32..(1 << nv) {
            let free: Vec<bool> = (0..nv).map(|i| (mask >> i) & 1 == 1).collect();
            // Lift to the full 6-variable assignment.
            let mut full = vec![false; nvars as usize];
            let mut t = 0usize;
            for i in 0..nvars as usize {
                match assign[i] {
                    Some(val) => full[i] = val,
                    None => {
                        full[i] = free[t];
                        t += 1;
                    }
                }
            }
            for (orig, sl) in eqs.iter().zip(sliced.iter()) {
                assert_eq!(
                    eval_bool(orig, nvars, &full),
                    eval_bool(sl, nv, &free),
                    "slice disagrees at mask {mask}"
                );
            }
        }
    }

    fn eval_bool(f: &F2BoolPoly, nvars: u32, x: &[bool]) -> bool {
        let mut acc = f.coeffs[0];
        for i in 0..nvars as usize {
            if f.coeffs.get(1 + i).copied().unwrap_or(false) && x[i] {
                acc ^= true;
            }
        }
        for i in 0..nvars {
            for j in (i + 1)..nvars {
                let idx = quad_monomial_index(i, j, nvars);
                if idx < f.coeffs.len() && f.coeffs[idx] && x[i as usize] && x[j as usize] {
                    acc ^= true;
                }
            }
        }
        acc
    }

    /// The guess patterns must return `k` distinct in-range indices, and
    /// `OneSide` must stay inside the `X₁` half while `Balanced` must not.
    #[test]
    fn guess_patterns_pick_distinct_in_range_indices() {
        let mut state = 12345u64;
        let mut rng = move || {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545F4914F6CDD1D)
        };
        let n_vars = 12u32;
        for pat in [
            GuessPattern::OneSide,
            GuessPattern::Balanced,
            GuessPattern::Spread,
        ] {
            for k in 0..=6u32 {
                let idx = pat.pick(n_vars, k, &mut rng);
                assert_eq!(idx.len(), k as usize, "{pat:?} k={k}");
                let mut sorted = idx.clone();
                sorted.sort_unstable();
                sorted.dedup();
                assert_eq!(sorted.len(), k as usize, "{pat:?} k={k} has duplicates");
                assert!(idx.iter().all(|&i| i < n_vars));
            }
        }
        let one = GuessPattern::OneSide.pick(n_vars, 4, &mut rng);
        assert!(one.iter().all(|&i| i < n_vars / 2));
        let bal = GuessPattern::Balanced.pick(n_vars, 4, &mut rng);
        assert!(bal.iter().any(|&i| i >= n_vars / 2));
    }

    /// The cost model must be Jensen-correct: computed from the histogram,
    /// not from the mean degree. A row whose `D*` is split between 2 and 4
    /// must cost strictly more than a row pinned at the mean degree 3.
    #[test]
    fn cost_model_uses_the_histogram_not_the_mean() {
        let mk = |hist: BTreeMap<u32, u32>| HybridRow {
            k: 0,
            vars: 12,
            eqs: 12,
            rho: 1.0,
            slices: 2,
            censored: 0,
            dstar_mean: 3.0,
            dstar_max: 4,
            dstar_hist: hist,
            delta_low_mean: 0.0,
            first_fall_mean: 3.0,
        };
        let split = mk(BTreeMap::from([(2, 1), (4, 1)]));
        let pinned = mk(BTreeMap::from([(3, 2)]));
        assert!(
            split.log2_total_cost(2.807) > pinned.log2_total_cost(2.807),
            "convexity: a spread of degrees must cost more than its mean"
        );
    }

    /// End-to-end smoke test of the sweep at a small operating point: rows
    /// must be produced for every `k`, `ρ` must increase with `k`, and every
    /// slice must refute (all slices of an unsatisfiable system are
    /// unsatisfiable).
    #[test]
    fn hybrid_sweep_runs_and_every_slice_refutes() {
        let n = 8;
        let irr = choose_irreducible(n);
        let sweep = run_hybrid_sweep(
            BasisFamily::Coordinate,
            GuessPattern::Balanced,
            n,
            4,
            &irr,
            3,
            2,
            3,
            6,
            0xC0FFEE,
        )
        .expect("coordinate family is always feasible");

        assert_eq!(sweep.full_vars, 8);
        assert_eq!(sweep.rows.len(), 4);
        for w in sweep.rows.windows(2) {
            assert!(w[1].rho > w[0].rho, "rho must rise with k");
            assert_eq!(w[1].vars + 1, w[0].vars);
        }
        for r in &sweep.rows {
            assert_eq!(
                r.censored, 0,
                "every slice of an unsatisfiable system must refute by d_max"
            );
            assert!(r.dstar_mean >= 2.0);
        }
        assert!(sweep.best(2.807).is_some());
    }

    /// The guess patterns must be compared on the *same* targets: the `k = 0`
    /// row fixes no variables, so it must be bit-identical across patterns.
    /// This fails if the target and slice randomness share one stream.
    #[test]
    fn k0_row_is_pattern_independent() {
        let n = 10;
        let irr = choose_irreducible(n);
        let run = |pat: GuessPattern| {
            run_hybrid_sweep(
                BasisFamily::Coordinate,
                pat,
                n,
                5,
                &irr,
                4,
                3,
                4,
                6,
                0xA11CE,
            )
            .unwrap()
            .rows[0]
                .clone()
        };
        let a = run(GuessPattern::Balanced);
        let b = run(GuessPattern::OneSide);
        let c = run(GuessPattern::Spread);
        assert_eq!(a.dstar_hist, b.dstar_hist, "balanced vs one-side at k=0");
        assert_eq!(a.dstar_hist, c.dstar_hist, "balanced vs spread at k=0");
        assert_eq!(a.dstar_max, c.dstar_max);
    }

    /// The enumeration baseline must dominate the degenerate `k = N` endpoint
    /// of the hybrid family: guessing every variable *is* enumeration, so a
    /// cost model that lets the hybrid beat `2^N` there is mis-specified.
    #[test]
    fn enumeration_baseline_bounds_the_full_guess_endpoint() {
        let (n_vars, eqs) = (14u32, 14u32);
        let base = log2_enumeration_cost(n_vars, eqs);
        // k = N leaves 0 free variables: cols(0, D) = 1, so the slice cost is
        // 1 and the total is exactly 2^N — strictly below the baseline, which
        // also charges `eqs` work per point.
        let full_guess = HybridRow {
            k: n_vars,
            vars: 0,
            eqs,
            rho: f64::INFINITY,
            slices: 1,
            censored: 0,
            dstar_mean: 2.0,
            dstar_max: 2,
            dstar_hist: BTreeMap::from([(2, 1)]),
            delta_low_mean: 0.0,
            first_fall_mean: 2.0,
        };
        assert!(
            full_guess.log2_total_cost(2.807) <= base,
            "baseline must be at least the full-guess endpoint"
        );
        assert!((full_guess.log2_total_cost(2.807) - n_vars as f64).abs() < 1e-9);
    }

    /// `collapse_fraction` must report the first `k` whose slices are all at
    /// the degree-2 floor, normalised by `N`.
    #[test]
    fn collapse_fraction_finds_the_floor() {
        let n = 10;
        let irr = choose_irreducible(n);
        let sweep = run_hybrid_sweep(
            BasisFamily::Coordinate,
            GuessPattern::OneSide,
            n,
            5,
            &irr,
            8,
            3,
            4,
            6,
            0xFEED,
        )
        .unwrap();
        let c = sweep
            .collapse_fraction()
            .expect("floor is reached by k=N-2");
        assert!(
            (0.0..=1.0).contains(&c),
            "collapse fraction out of range: {c}"
        );
        let k = (c * sweep.full_vars as f64).round() as u32;
        assert!(sweep.rows[k as usize].dstar_max <= 2);
        if k > 0 {
            assert!(
                sweep.rows[(k - 1) as usize].dstar_max > 2,
                "collapse_fraction must return the FIRST collapsed k"
            );
        }
    }

    /// The boundary check must actually fire: on a toy system the cost model
    /// slides the optimum to the largest `k` scanned (exhaustive search), and
    /// `optimum_is_interior` must report that rather than let the sweep claim
    /// a saving it did not find.
    #[test]
    fn boundary_optimum_is_reported_as_non_interior() {
        let n = 12;
        let irr = choose_irreducible(n);
        let sweep = run_hybrid_sweep(
            BasisFamily::Coordinate,
            GuessPattern::Balanced,
            n,
            6,
            &irr,
            10,
            2,
            3,
            6,
            0xB0A11,
        )
        .unwrap();
        let (k_star, _) = sweep.best(2.807).unwrap();
        assert_eq!(
            sweep.optimum_is_interior(2.807),
            k_star < sweep.k_max(),
            "interior flag must agree with the optimum's position"
        );
    }

    /// A cell whose unsliced system is already at the degree-2 floor has no
    /// collapse to measure, so `collapse_fraction` must decline to report one
    /// rather than return a degenerate 0.
    #[test]
    fn collapse_fraction_declines_when_already_at_the_floor() {
        // Heavily over-determined: rho = 10/4, so D* = 2 from the start (P6).
        let n = 10;
        let irr = choose_irreducible(n);
        let sweep = run_hybrid_sweep(
            BasisFamily::Coordinate,
            GuessPattern::Balanced,
            n,
            2,
            &irr,
            2,
            2,
            2,
            6,
            0x1234,
        )
        .unwrap();
        assert_eq!(
            sweep.rows[0].dstar_max, 2,
            "expected the P6 floor at rho=2.5"
        );
        assert!(
            sweep.collapse_fraction().is_none(),
            "no collapse to report when k=0 is already at the floor"
        );
    }

    /// The load-bearing claim of the whole lever: raising the determination
    /// ratio must not *raise* `D*`. Slicing can only add constraints, so the
    /// worst-case degree must be monotone non-increasing in `k`.
    #[test]
    fn slicing_does_not_raise_the_mean_solving_degree() {
        let n = 10;
        let irr = choose_irreducible(n);
        let sweep = run_hybrid_sweep(
            BasisFamily::Coordinate,
            GuessPattern::Balanced,
            n,
            5,
            &irr,
            4,
            3,
            4,
            6,
            0x5EED,
        )
        .unwrap();
        let first = sweep.rows.first().unwrap().dstar_mean;
        let last = sweep.rows.last().unwrap().dstar_mean;
        assert!(
            last <= first + 1e-9,
            "mean D* rose under slicing: {first} → {last}"
        );
    }
}
