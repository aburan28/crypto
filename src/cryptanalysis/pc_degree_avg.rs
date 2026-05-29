//! # Refutation-degree **distribution** over many non-decomposable targets.
//!
//! Averaging layer on top of [`crate::cryptanalysis::pc_degree_harness`].
//! A single `(n, n', b, x₃)` instance gives one noisy refutation degree
//! `D*`; this module samples **many** non-decomposable targets per
//! `(n, n')` and reports the `D*` *distribution* (min / mean / max plus a
//! histogram), which is what the proposal's prediction #1 actually talks
//! about (`RESEARCH_FFD_PROOF_COMPLEXITY.md` §4 item 1, §6).
//!
//! Only **non-decomposable** (genuinely unsatisfiable) targets are
//! counted — those are the instances for which a refutation exists and
//! `D*` is defined. Decomposable draws are skipped (and tallied as
//! `skipped_decomposable`) so the reported statistics are conditional on
//! "the target does not decompose over `V`", exactly the index-calculus
//! failure case whose cost the solving degree governs.
//!
//! ## Why a distribution, not a point
//!
//! At toy sizes the per-target `D*` jumps around (some targets admit a
//! degree-2 Nullstellensatz certificate, others need more), so a single
//! sample tells you little about the asymptotic trend. The **mean** and
//! **max** over a batch are the stable signals: prediction #1 says the
//! mean `D*` climbs with `n` while the first fall degree stays flat. This
//! module produces the numbers that test that, and the histogram exposes
//! the spread that a single sample hides.
//!
//! ## Cost
//!
//! Each target costs one `pc_degree_harness::measure_refutation_degree`
//! (a Macaulay rank per degree up to `D*`) plus one decomposability
//! brute-force over `|V|² = 2^{2n'}`. Comfortable for `n' ≤ 6` and
//! batches of a few hundred targets. For larger `n'` the decomposability
//! oracle dominates and should be replaced by an `S₃`-resultant solve.
//!
//! ## References
//!
//! See `pc_degree_harness` and `RESEARCH_FFD_PROOF_COMPLEXITY.md`.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::ffd_harness::{choose_irreducible, random_nonzero_f2m};
use crate::cryptanalysis::pc_degree_harness::{is_decomposable, measure_refutation_degree};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::collections::BTreeMap;

// ── Public API ──────────────────────────────────────────────────────

/// Summary statistics of an integer-valued degree over a batch of
/// samples: the count, none-count, min/max, mean, and a histogram
/// (degree → frequency).
#[derive(Clone, Debug, Default)]
pub struct DegreeStats {
    /// Number of samples with a *defined* degree (finite, not `None`).
    pub n_defined: u32,
    /// Number of samples where the degree was `None` (no event observed
    /// within the degree budget).
    pub n_none: u32,
    /// Minimum observed degree (over defined samples).
    pub min: Option<u32>,
    /// Maximum observed degree (over defined samples).
    pub max: Option<u32>,
    /// Mean of the defined degrees.
    pub mean: Option<f64>,
    /// Histogram: degree → number of samples with that degree.
    pub histogram: BTreeMap<u32, u32>,
}

impl DegreeStats {
    fn accumulate(&mut self, value: Option<u32>) {
        match value {
            Some(d) => {
                self.n_defined += 1;
                self.min = Some(self.min.map_or(d, |m| m.min(d)));
                self.max = Some(self.max.map_or(d, |m| m.max(d)));
                *self.histogram.entry(d).or_insert(0) += 1;
            }
            None => self.n_none += 1,
        }
    }

    fn finalize(&mut self) {
        if self.n_defined > 0 {
            let sum: u64 = self
                .histogram
                .iter()
                .map(|(&d, &c)| d as u64 * c as u64)
                .sum();
            self.mean = Some(sum as f64 / self.n_defined as f64);
        }
    }

    /// Compact histogram rendering, e.g. `2:3 3:11 4:2`.
    pub fn histogram_str(&self) -> String {
        if self.histogram.is_empty() {
            return "—".into();
        }
        self.histogram
            .iter()
            .map(|(d, c)| format!("{d}:{c}"))
            .collect::<Vec<_>>()
            .join(" ")
    }
}

/// One row of the averaged sweep: a `(n, n')` cell with the `D*` and
/// first-fall distributions over a batch of non-decomposable targets.
#[derive(Clone, Debug)]
pub struct AvgRow {
    pub n: u32,
    pub n_sub: u32,
    pub num_vars: u32,
    pub num_eqs: u32,
    /// Non-decomposable targets actually measured.
    pub n_targets: u32,
    /// Decomposable draws skipped while hunting for non-decomposable ones.
    pub skipped_decomposable: u32,
    /// Distribution of the refutation degree `D*`.
    pub refutation: DegreeStats,
    /// Distribution of the first-fall degree, for the same instances.
    pub first_fall: DegreeStats,
}

/// Measure the `D*` distribution over up to `targets` non-decomposable
/// targets for a fixed `(n, n', b)`. Draws random `x₃` with `rng`,
/// skipping decomposable ones, until `targets` non-decomposable instances
/// are measured or `max_attempts` draws are exhausted.
pub fn measure_avg(
    n: u32,
    n_sub: u32,
    irr: &IrreduciblePoly,
    b: &F2mElement,
    d_max: u32,
    targets: u32,
    max_attempts: u32,
    rng: &mut StdRng,
) -> AvgRow {
    let mut refutation = DegreeStats::default();
    let mut first_fall = DegreeStats::default();
    let mut n_targets = 0u32;
    let mut skipped = 0u32;
    let mut attempts = 0u32;

    while n_targets < targets && attempts < max_attempts {
        attempts += 1;
        let x3 = random_nonzero_f2m(rng, n);
        if is_decomposable(n, n_sub, irr, b, &x3) {
            skipped += 1;
            continue;
        }
        let row = measure_refutation_degree(n, n_sub, irr, b, &x3, d_max);
        debug_assert!(!row.decomposable);
        refutation.accumulate(row.refutation_degree);
        first_fall.accumulate(row.first_fall);
        n_targets += 1;
    }

    refutation.finalize();
    first_fall.finalize();

    AvgRow {
        n,
        n_sub,
        num_vars: 2 * n_sub,
        num_eqs: n,
        n_targets,
        skipped_decomposable: skipped,
        refutation,
        first_fall,
    }
}

/// Run the averaged sweep over `n_range`, with `n' = ⌊n/2⌋`, `targets`
/// non-decomposable instances per `n`, Macaulay degree cap `d_max`, and a
/// reproducible `seed`. `max_attempts` bounds the random search per `n`
/// (so small fields that lack enough non-decomposable targets terminate).
pub fn run_pc_avg_sweep(
    n_range: std::ops::RangeInclusive<u32>,
    d_max: u32,
    targets: u32,
    seed: u64,
) -> Vec<AvgRow> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out = Vec::new();
    for n in n_range {
        let n_sub = (n / 2).max(1);
        let irr = choose_irreducible(n);
        let b = random_nonzero_f2m(&mut rng, n);
        // Allow generous head-room over `targets` so decomposable draws
        // don't starve the batch.
        let max_attempts = targets.saturating_mul(20).max(64);
        out.push(measure_avg(
            n,
            n_sub,
            &irr,
            &b,
            d_max,
            targets,
            max_attempts,
            &mut rng,
        ));
    }
    out
}

// ── Pretty printing ─────────────────────────────────────────────────

pub fn print_pc_avg_sweep(rows: &[AvgRow]) {
    println!();
    println!(
        "{:>3} {:>3} {:>5} {:>4} {:>5} | {:>13} | {:>17} | {:>4}",
        "n", "n'", "vars", "eqs", "N", "first-fall μ", "D* min/μ/max", "none"
    );
    println!("{}", "─".repeat(72));
    for r in rows {
        let ff_mean = r
            .first_fall
            .mean
            .map(|m| format!("{m:.2}"))
            .unwrap_or_else(|| "—".into());
        let ds = &r.refutation;
        let dstar = match (ds.min, ds.mean, ds.max) {
            (Some(lo), Some(mu), Some(hi)) => format!("{lo}/{mu:.2}/{hi}"),
            _ => "—".into(),
        };
        println!(
            "{:>3} {:>3} {:>5} {:>4} {:>5} | {:>13} | {:>17} | {:>4}",
            r.n,
            r.n_sub,
            r.num_vars,
            r.num_eqs,
            r.n_targets,
            ff_mean,
            dstar,
            ds.n_none,
        );
    }
    println!();
    println!("D* histograms (degree:count):");
    for r in rows {
        println!(
            "  n={:>2} n'={:<2} → {}   [skipped {} decomposable draws]",
            r.n,
            r.n_sub,
            r.refutation.histogram_str(),
            r.skipped_decomposable,
        );
    }
    println!();
    println!(
        "Reading: `D* min/μ/max` is the refutation-degree distribution over\n\
         N non-decomposable (unsatisfiable) targets — the stable form of the\n\
         Polynomial-Calculus / last-fall degree. Prediction #1: mean D* climbs\n\
         with n while the first-fall mean stays roughly flat. `none` counts\n\
         instances that did not refute within the degree budget (raise d_max\n\
         if non-zero)."
    );
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// `DegreeStats` computes min/max/mean/histogram correctly and counts
    /// `None`s separately.
    #[test]
    fn degree_stats_arithmetic() {
        let mut s = DegreeStats::default();
        for v in [Some(2u32), Some(4), Some(2), None, Some(3)] {
            s.accumulate(v);
        }
        s.finalize();
        assert_eq!(s.n_defined, 4);
        assert_eq!(s.n_none, 1);
        assert_eq!(s.min, Some(2));
        assert_eq!(s.max, Some(4));
        assert_eq!(s.mean, Some((2.0 + 4.0 + 2.0 + 3.0) / 4.0));
        assert_eq!(s.histogram.get(&2), Some(&2));
        assert_eq!(s.histogram.get(&3), Some(&1));
        assert_eq!(s.histogram.get(&4), Some(&1));
        assert_eq!(s.histogram_str(), "2:2 3:1 4:1");
    }

    /// Empty stats stay empty and render as `—`.
    #[test]
    fn degree_stats_empty() {
        let mut s = DegreeStats::default();
        s.finalize();
        assert_eq!(s.n_defined, 0);
        assert!(s.mean.is_none());
        assert_eq!(s.histogram_str(), "—");
    }

    /// **Every measured target in a batch is non-decomposable**, so the
    /// reported `D*` distribution is over genuinely unsatisfiable
    /// instances, and the histogram sample-count matches `n_targets`.
    #[test]
    fn batch_only_counts_nondecomposable_and_is_consistent() {
        let n = 4;
        let n_sub = 2;
        let irr = choose_irreducible(n);
        let mut rng = StdRng::seed_from_u64(0x_AA_BB);
        let b = random_nonzero_f2m(&mut rng, n);
        let row = measure_avg(n, n_sub, &irr, &b, 2 * n_sub + 1, 12, 4000, &mut rng);
        assert_eq!(row.num_eqs, n);
        assert_eq!(row.num_vars, 2 * n_sub);
        // Histogram totals (defined) + none-count == targets measured.
        let hist_total: u32 = row.refutation.histogram.values().sum();
        assert_eq!(hist_total, row.refutation.n_defined);
        assert_eq!(
            row.refutation.n_defined + row.refutation.n_none,
            row.n_targets
        );
        // At n=4/n'=2 the field is small but at least one non-decomposable
        // target should be found.
        assert!(row.n_targets >= 1, "no non-decomposable target found");
        // Non-decomposable ⇒ a refutation exists, so within a generous
        // degree budget there should be no `none`s.
        assert_eq!(row.refutation.n_none, 0, "non-decomposable target failed to refute");
    }

    /// The averaged sweep runs end-to-end and produces one row per `n`.
    #[test]
    fn avg_sweep_runs() {
        let rows = run_pc_avg_sweep(4..=6, 9, 8, 0x_A6_DE);
        assert_eq!(rows.len(), 3);
        for r in &rows {
            assert_eq!(r.num_eqs, r.n);
            // Mean is defined whenever at least one target refuted.
            if r.refutation.n_defined > 0 {
                assert!(r.refutation.mean.is_some());
                assert!(r.refutation.min.unwrap() <= r.refutation.max.unwrap());
            }
        }
    }
}
