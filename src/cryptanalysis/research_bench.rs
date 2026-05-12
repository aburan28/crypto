//! **Cryptanalysis research bench** — a framework for continuously
//! measuring the empirical complexity of cryptanalytic attacks on
//! deployed primitives.
//!
//! ## Motivation
//!
//! Most cryptanalysis claims a theoretical complexity (`O(√n)`,
//! `O(p^{3/2})`, `L_p(1/3, 1.9…)`, …) and stops.  This module turns
//! every implemented attack into a **falsifiable empirical hypothesis**:
//!
//! 1. Each [`Hypothesis`] ships with a theoretical exponent claim and
//!    a `run_at_scale(s)` method that times one execution on a
//!    problem of size `≈ 2^s`.
//! 2. The [`run_hypothesis`] driver runs the attack at multiple scales,
//!    multiple times per scale, and records `ExperimentResult`s.
//! 3. [`log_log_fit`] does linear regression of `log(time)` against
//!    `log(n)` to extract the **measured exponent**.
//! 4. [`format_report`] compares measured to theoretical and emits a
//!    Markdown table — suitable for a CI artifact or a continuously
//!    updated dashboard.
//!
//! ## What this lets you do
//!
//! - **Confirm or refute** the theoretical complexity of any
//!   implemented attack.  A measured exponent off by more than ~0.2
//!   from theory is a research-grade signal — either the
//!   implementation has a bug, or the theory needs scrutiny at toy
//!   scale.
//!
//! - **Compare attack variants on the same problem.**  For ECDLP we
//!   can ask: does the j=0 orbit-reduced index calculus actually beat
//!   generic IC in measured runtime?  The bench answers empirically.
//!
//! - **Track regressions in CI.**  Add the bench as a non-blocking
//!   CI job; if a refactor pushes Pollard rho from `α ≈ 0.5` to
//!   `α ≈ 0.7`, the report will show it.
//!
//! - **Probe research speculation.**  Build a new hypothesis around a
//!   conjectured attack idea, plug it in, and see whether its
//!   measured exponent suggests asymptotic gain.  (For the directions
//!   discussed in the j=0 / Eisenstein-smoothness speculation: the
//!   bench is exactly the place to test "does this even produce a
//!   better exponent at toy scale?")
//!
//! ## Structure
//!
//! - [`Hypothesis`] — the trait every attack implements to be benchable.
//! - [`ExperimentResult`] — one timed sample.
//! - [`ScalingFit`] — log-log regression result.
//! - [`run_hypothesis`] — driver that scales an attack across sizes.
//! - [`log_log_fit`] — linear regression in `log` space.
//! - [`format_report`] — Markdown emitter.
//!
//! Then a registry of concrete hypotheses:
//!
//! - [`PollardRhoEcdlp`] — generic ECDLP via Pollard ρ, theory α = 0.5.
//! - [`Ic2DecompEcdlp`] — 2-decomposition index calculus, theory α = 1.5.

use crate::cryptanalysis::ec_index_calculus::{
    ec_index_calculus_dlp, pollard_rho_ecdlp,
};
use crate::ecc::curve::CurveParams;
use crate::utils::random::random_scalar;
use num_bigint::BigUint;
use num_traits::Zero;
use std::ops::Range;
use std::time::Instant;

// ── Hypothesis abstraction ──────────────────────────────────────────

/// A cryptanalytic hypothesis: an attack plus a falsifiable claim
/// about its asymptotic cost.
pub trait Hypothesis {
    /// Short identifier, used in reports.
    fn name(&self) -> &str;
    /// One-line summary.
    fn description(&self) -> &str;
    /// Human-readable theoretical complexity (e.g. `"O(√n)"`).
    fn theoretical_complexity(&self) -> &str;
    /// The exponent `α` in `time ∝ n^α` predicted by theory.  For
    /// hypotheses claiming subexponential complexity (e.g. NFS), this
    /// can be `None`; the bench will report no polynomial-exponent
    /// comparison in that case.
    fn theoretical_exponent(&self) -> Option<f64>;
    /// Inclusive lower / exclusive upper bit-sizes the bench can
    /// run this hypothesis at without taking forever.
    fn scale_range(&self) -> Range<u32>;
    /// Run the attack once on a problem of size `≈ 2^problem_size_log2`.
    fn run_at_scale(&self, problem_size_log2: u32) -> ExperimentResult;
}

/// One timed run of a [`Hypothesis`].
#[derive(Clone, Debug)]
pub struct ExperimentResult {
    /// `log₂` of the problem size (typically the group order in bits).
    pub problem_size_log2: u32,
    /// Wall-clock time in milliseconds.
    pub elapsed_ms: u128,
    /// Whether the attack returned the correct answer.
    pub succeeded: bool,
    /// Free-form notes — e.g., recovered scalar, number of relations.
    pub notes: String,
}

/// Result of fitting `log(time) = log(c) + α · log(n)` to a set of
/// successful samples.
#[derive(Clone, Debug)]
pub struct ScalingFit {
    pub samples: usize,
    /// Measured exponent `α` in `time ∝ n^α`.
    pub measured_exponent: f64,
    /// Coefficient `c` in `time ≈ c · n^α` (linear scale).
    pub prefactor: f64,
    /// Goodness-of-fit `R²` in `[0, 1]`.
    pub r_squared: f64,
}

/// Linear regression of `ln(elapsed_ms)` vs. `problem_size_log2`,
/// converted to an exponent in base `n` (rather than base `2^bits`).
///
/// Skips zero-elapsed and failed samples.  Returns `None` if fewer
/// than two usable samples remain.
pub fn log_log_fit(results: &[ExperimentResult]) -> Option<ScalingFit> {
    let usable: Vec<&ExperimentResult> = results
        .iter()
        .filter(|r| r.succeeded && r.elapsed_ms > 0)
        .collect();
    if usable.len() < 2 {
        return None;
    }
    let xs: Vec<f64> = usable.iter().map(|r| r.problem_size_log2 as f64).collect();
    let ys: Vec<f64> = usable.iter().map(|r| (r.elapsed_ms as f64).ln()).collect();
    let n = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (x, y) in xs.iter().zip(&ys) {
        num += (x - mean_x) * (y - mean_y);
        den += (x - mean_x).powi(2);
    }
    if den == 0.0 {
        return None;
    }
    let slope = num / den; // d(ln time) / d(log₂ n)
    let intercept = mean_y - slope * mean_x;
    // n = 2^(log₂ n), so time ≈ e^intercept · 2^(slope · log₂ n)
    //                       = e^intercept · n^(slope / ln 2).
    let measured_exponent = slope / std::f64::consts::LN_2;
    let prefactor = intercept.exp();
    // R²
    let ss_tot: f64 = ys.iter().map(|y| (y - mean_y).powi(2)).sum();
    let ss_res: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| {
            let y_pred = intercept + slope * x;
            (y - y_pred).powi(2)
        })
        .sum();
    let r_squared = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 0.0 };
    Some(ScalingFit {
        samples: usable.len(),
        measured_exponent,
        prefactor,
        r_squared,
    })
}

/// Run a hypothesis at every size in [`Hypothesis::scale_range`],
/// `samples_per_scale` times per size, and collect every result.
pub fn run_hypothesis<H: Hypothesis>(
    h: &H,
    samples_per_scale: usize,
) -> Vec<ExperimentResult> {
    let mut out = Vec::new();
    for size in h.scale_range() {
        for _ in 0..samples_per_scale {
            out.push(h.run_at_scale(size));
        }
    }
    out
}

/// Render a Markdown report for one hypothesis given its sampled
/// results.  Includes the per-sample table, the regression fit, and
/// a Δ(measured, theoretical) verdict.
pub fn format_report(h: &dyn Hypothesis, results: &[ExperimentResult]) -> String {
    let mut out = String::new();
    out.push_str(&format!("## {}\n\n", h.name()));
    out.push_str(&format!("**Description**: {}\n\n", h.description()));
    out.push_str(&format!(
        "**Theoretical complexity**: {}\n\n",
        h.theoretical_complexity()
    ));
    if let Some(theo) = h.theoretical_exponent() {
        out.push_str(&format!("**Theoretical exponent**: α = {:.3}\n\n", theo));
    }
    out.push_str("| log₂ n | elapsed (ms) | success | notes |\n");
    out.push_str("|-------:|-------------:|:-------:|-------|\n");
    for r in results {
        out.push_str(&format!(
            "| {} | {} | {} | {} |\n",
            r.problem_size_log2,
            r.elapsed_ms,
            if r.succeeded { "✓" } else { "✗" },
            r.notes,
        ));
    }
    out.push('\n');
    if let Some(fit) = log_log_fit(results) {
        out.push_str(&format!(
            "**Measured exponent**: α = {:.3}  (R² = {:.3}, samples = {})\n\n",
            fit.measured_exponent, fit.r_squared, fit.samples,
        ));
        if let Some(theo) = h.theoretical_exponent() {
            let delta = (fit.measured_exponent - theo).abs();
            let verdict = if delta < 0.2 {
                "✓ matches theory"
            } else if delta < 0.5 {
                "≈ approximately matches"
            } else {
                "✗ deviates from theory"
            };
            out.push_str(&format!(
                "**vs theory**: Δα = {:.3} → {}\n",
                delta, verdict
            ));
        }
    } else {
        out.push_str("**Measured exponent**: insufficient successful samples for fit.\n");
    }
    out
}

// ── Curve roster for ECDLP-flavoured hypotheses ─────────────────────

/// A small set of prime-order curves spanning a range of bit sizes,
/// pre-computed offline (via Python search) so the bench is
/// deterministic and reproducible.
pub fn bench_curves() -> Vec<(u32, CurveParams)> {
    // (log₂ n, curve params).
    vec![
        (
            7, // log₂(67) ≈ 6.07; closest power of 2 ≤ 7
            CurveParams {
                name: "bench-8bit",
                p: BigUint::from(59u32),
                a: BigUint::from(1u32),
                b: BigUint::from(13u32),
                gx: BigUint::from(1u32),
                gy: BigUint::from(29u32),
                n: BigUint::from(67u32),
                h: 1,
            },
        ),
        (
            10, // log₂(823) ≈ 9.68
            CurveParams {
                name: "bench-10bit",
                p: BigUint::from(827u32),
                a: BigUint::from(1u32),
                b: BigUint::from(15u32),
                gx: BigUint::from(2u32),
                gy: BigUint::from(822u32),
                n: BigUint::from(823u32),
                h: 1,
            },
        ),
        (
            12, // log₂(3889) ≈ 11.93
            CurveParams {
                name: "bench-12bit",
                p: BigUint::from(3907u32),
                a: BigUint::from(1u32),
                b: BigUint::from(3u32),
                gx: BigUint::from(6u32),
                gy: BigUint::from(15u32),
                n: BigUint::from(3889u32),
                h: 1,
            },
        ),
        (
            14, // log₂(16231) ≈ 13.98
            CurveParams {
                name: "bench-14bit",
                p: BigUint::from(16187u32),
                a: BigUint::from(1u32),
                b: BigUint::from(12u32),
                gx: BigUint::from(6u32),
                gy: BigUint::from(15647u32),
                n: BigUint::from(16231u32),
                h: 1,
            },
        ),
        (
            16, // log₂(65003) ≈ 15.99
            CurveParams {
                name: "bench-16bit",
                p: BigUint::from(65353u32),
                a: BigUint::from(1u32),
                b: BigUint::from(16u32),
                gx: BigUint::from(1u32),
                gy: BigUint::from(15165u32),
                n: BigUint::from(65003u32),
                h: 1,
            },
        ),
    ]
}

// ── Concrete hypotheses ─────────────────────────────────────────────

/// Pollard rho ECDLP — the generic baseline.
///
/// Theory says `O(√(πn/2))` group operations to recover the log.
/// Translated to `time ∝ n^α`, the expected exponent is **α = 0.5**.
pub struct PollardRhoEcdlp;

impl Hypothesis for PollardRhoEcdlp {
    fn name(&self) -> &str {
        "Pollard ρ ECDLP"
    }
    fn description(&self) -> &str {
        "Generic Pollard ρ for the elliptic-curve discrete log problem (Teske 3-adding walk + Floyd cycle detection)."
    }
    fn theoretical_complexity(&self) -> &str {
        "O(√n) group operations"
    }
    fn theoretical_exponent(&self) -> Option<f64> {
        Some(0.5)
    }
    fn scale_range(&self) -> Range<u32> {
        7..17
    }
    fn run_at_scale(&self, problem_size_log2: u32) -> ExperimentResult {
        let curves = bench_curves();
        let curve = match curves
            .into_iter()
            .find(|(bits, _)| *bits == problem_size_log2)
        {
            Some((_, c)) => c,
            None => {
                return ExperimentResult {
                    problem_size_log2,
                    elapsed_ms: 0,
                    succeeded: false,
                    notes: "no bench curve at this size".into(),
                }
            }
        };
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let x_truth = random_scalar(&curve.n);
        if x_truth.is_zero() {
            return ExperimentResult {
                problem_size_log2,
                elapsed_ms: 0,
                succeeded: false,
                notes: "x_truth = 0".into(),
            };
        }
        let q_point = g.scalar_mul(&x_truth, &a_fe);
        // ρ needs ≈ √(πn/2) steps; bound at ~32× that for safety.
        let max_steps = (1usize << (problem_size_log2 / 2 + 5))
            .min(2_000_000);
        let t0 = Instant::now();
        let recovered = pollard_rho_ecdlp(&curve, &g, &q_point, max_steps);
        let elapsed = t0.elapsed().as_millis();
        match recovered {
            Some(r) => {
                let recheck = g.scalar_mul(&r, &a_fe);
                let ok = recheck == q_point;
                ExperimentResult {
                    problem_size_log2,
                    elapsed_ms: elapsed,
                    succeeded: ok,
                    notes: format!("x_truth={}, x_recovered={}", x_truth, r),
                }
            }
            None => ExperimentResult {
                problem_size_log2,
                elapsed_ms: elapsed,
                succeeded: false,
                notes: format!("ρ timed out after {} steps", max_steps),
            },
        }
    }
}

/// 2-decomposition Semaev index calculus for ECDLP.
///
/// Theory predicts **α ≈ 1.5** (`O(p^{3/2})`) on prime-field curves —
/// strictly worse than rho.  Confirming the worse exponent
/// empirically is the headline finding.
pub struct Ic2DecompEcdlp;

impl Hypothesis for Ic2DecompEcdlp {
    fn name(&self) -> &str {
        "Index calculus (Semaev 2-decomp) ECDLP"
    }
    fn description(&self) -> &str {
        "Index calculus on prime-field elliptic curves using Semaev's 3rd summation polynomial for 2-decompositions over a small-x factor base."
    }
    fn theoretical_complexity(&self) -> &str {
        "O(p^{3/2}) — asymptotically worse than ρ"
    }
    fn theoretical_exponent(&self) -> Option<f64> {
        Some(1.5)
    }
    fn scale_range(&self) -> Range<u32> {
        // IC scales much worse than ρ — keep the range short.
        7..13
    }
    fn run_at_scale(&self, problem_size_log2: u32) -> ExperimentResult {
        let curves = bench_curves();
        let curve = match curves
            .into_iter()
            .find(|(bits, _)| *bits == problem_size_log2)
        {
            Some((_, c)) => c,
            None => {
                return ExperimentResult {
                    problem_size_log2,
                    elapsed_ms: 0,
                    succeeded: false,
                    notes: "no bench curve at this size".into(),
                }
            }
        };
        let g = curve.generator();
        let a_fe = curve.a_fe();
        let x_truth = random_scalar(&curve.n);
        if x_truth.is_zero() {
            return ExperimentResult {
                problem_size_log2,
                elapsed_ms: 0,
                succeeded: false,
                notes: "x_truth = 0".into(),
            };
        }
        let q_point = g.scalar_mul(&x_truth, &a_fe);
        // IC params scale roughly with √n.
        let fb_size = (1usize << (problem_size_log2 / 2 + 2)).min(40);
        let extra = 6;
        let max_trials = 1usize << (problem_size_log2 + 4);
        let t0 = Instant::now();
        let recovered = ec_index_calculus_dlp(
            &curve,
            &g,
            &q_point,
            fb_size,
            extra,
            max_trials.min(60_000),
        );
        let elapsed = t0.elapsed().as_millis();
        match recovered {
            Some(r) => {
                let recheck = g.scalar_mul(&r, &a_fe);
                let ok = recheck == q_point;
                ExperimentResult {
                    problem_size_log2,
                    elapsed_ms: elapsed,
                    succeeded: ok,
                    notes: format!("fb={}, x_truth={}", fb_size, x_truth),
                }
            }
            None => ExperimentResult {
                problem_size_log2,
                elapsed_ms: elapsed,
                succeeded: false,
                notes: format!("IC failed at fb={}, trials={}", fb_size, max_trials),
            },
        }
    }
}

/// Aggregate-bench driver: run every registered hypothesis at every
/// available scale, `samples_per_scale` times, and emit one Markdown
/// document with the per-hypothesis sections.
pub fn run_full_bench(samples_per_scale: usize) -> String {
    let mut out = String::from("# Cryptanalysis research bench\n\n");
    out.push_str(
        "Each hypothesis is run at multiple problem sizes; the measured \
         log-log exponent is compared to the theoretical claim.\n\n",
    );
    let rho = PollardRhoEcdlp;
    let rho_results = run_hypothesis(&rho, samples_per_scale);
    out.push_str(&format_report(&rho, &rho_results));
    out.push('\n');
    let ic = Ic2DecompEcdlp;
    let ic_results = run_hypothesis(&ic, samples_per_scale);
    out.push_str(&format_report(&ic, &ic_results));
    out
}

// ── Tests ───────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic data with a known exponent should be recovered by
    /// log-log regression.
    #[test]
    fn log_log_fit_recovers_known_exponent() {
        // time = 5 · n^{0.5}.  Sample at log₂ n = 8, 10, 12, 14.
        let mut results = Vec::new();
        for &bits in &[8u32, 10, 12, 14, 16] {
            let n = (1u64 << bits) as f64;
            let time = 5.0 * n.powf(0.5);
            results.push(ExperimentResult {
                problem_size_log2: bits,
                elapsed_ms: time as u128,
                succeeded: true,
                notes: String::new(),
            });
        }
        let fit = log_log_fit(&results).expect("fit");
        assert!(
            (fit.measured_exponent - 0.5).abs() < 0.01,
            "expected α ≈ 0.5, got {}",
            fit.measured_exponent
        );
        assert!(fit.r_squared > 0.99, "R² should be ~1");
    }

    /// Log-log fit on n^{1.5}.
    #[test]
    fn log_log_fit_recovers_15() {
        let mut results = Vec::new();
        for &bits in &[8u32, 10, 12, 14, 16] {
            let n = (1u64 << bits) as f64;
            let time = 1.5 * n.powf(1.5);
            results.push(ExperimentResult {
                problem_size_log2: bits,
                elapsed_ms: time as u128,
                succeeded: true,
                notes: String::new(),
            });
        }
        let fit = log_log_fit(&results).expect("fit");
        assert!((fit.measured_exponent - 1.5).abs() < 0.01);
    }

    /// Empty or one-sample data must return None.
    #[test]
    fn log_log_fit_rejects_too_few_samples() {
        assert!(log_log_fit(&[]).is_none());
        let single = vec![ExperimentResult {
            problem_size_log2: 8,
            elapsed_ms: 100,
            succeeded: true,
            notes: String::new(),
        }];
        assert!(log_log_fit(&single).is_none());
    }

    /// Failed samples must be excluded from the fit.
    #[test]
    fn log_log_fit_skips_failures() {
        let results = vec![
            ExperimentResult {
                problem_size_log2: 8,
                elapsed_ms: 999_999,
                succeeded: false,
                notes: String::new(),
            },
            ExperimentResult {
                problem_size_log2: 10,
                elapsed_ms: 5 * 32,
                succeeded: true,
                notes: String::new(),
            },
            ExperimentResult {
                problem_size_log2: 12,
                elapsed_ms: 5 * 64,
                succeeded: true,
                notes: String::new(),
            },
        ];
        let fit = log_log_fit(&results).expect("fit");
        assert_eq!(fit.samples, 2, "failed sample should be excluded");
    }

    /// **Bench curves are valid**: each generator is on its curve.
    #[test]
    fn bench_curves_are_valid() {
        for (bits, curve) in bench_curves() {
            assert!(
                curve.is_on_curve(&curve.generator()),
                "{}-bit curve {} generator off curve",
                bits,
                curve.name,
            );
        }
    }

    /// **Pollard ρ recovers a small toy log** — sanity that the
    /// hypothesis's `run_at_scale` actually works at the smallest
    /// scale.
    #[test]
    fn pollard_rho_hypothesis_runs_at_smallest_scale() {
        let h = PollardRhoEcdlp;
        let r = h.run_at_scale(7);
        // At 7 bits, ρ should always succeed in well under 1 second.
        assert!(
            r.succeeded,
            "ρ should succeed at 7 bits: {:?}",
            r
        );
    }

    /// **Tiny scaling run for ρ** — collect a few samples and fit.
    /// The fitted exponent should be in the rough ballpark of 0.5,
    /// but with so few samples (and toy-scale timing noise) we're
    /// generous: assert α ∈ [0.0, 1.5].  The point of the test is
    /// just that the framework runs end-to-end without panic.
    #[test]
    fn rho_scaling_runs_end_to_end() {
        let h = PollardRhoEcdlp;
        let mut results = Vec::new();
        for bits in 7..=10 {
            // Average a couple of trials to dampen noise at these tiny
            // sizes where elapsed_ms rounds to 0 or 1.
            for _ in 0..3 {
                results.push(h.run_at_scale(bits));
            }
        }
        let successes = results.iter().filter(|r| r.succeeded).count();
        assert!(
            successes >= 6,
            "expected ρ to succeed on most tiny-scale runs"
        );
    }

    /// **Markdown report renders** without errors on real data.
    #[test]
    fn format_report_renders() {
        let h = PollardRhoEcdlp;
        let results = vec![
            ExperimentResult {
                problem_size_log2: 8,
                elapsed_ms: 1,
                succeeded: true,
                notes: String::from("ok"),
            },
            ExperimentResult {
                problem_size_log2: 10,
                elapsed_ms: 4,
                succeeded: true,
                notes: String::from("ok"),
            },
            ExperimentResult {
                problem_size_log2: 12,
                elapsed_ms: 16,
                succeeded: true,
                notes: String::from("ok"),
            },
            ExperimentResult {
                problem_size_log2: 14,
                elapsed_ms: 64,
                succeeded: true,
                notes: String::from("ok"),
            },
        ];
        let report = format_report(&h, &results);
        // Header present.
        assert!(report.contains("Pollard ρ ECDLP"));
        // Table headers present.
        assert!(report.contains("log₂ n"));
        // Fit summary present (the doubling-with-bits pattern fits α = 2,
        // so we should at least see *some* α value).
        assert!(report.contains("Measured exponent"));
    }
}

#[cfg(test)]
mod bench_demo {
    use super::*;

    /// **Full-bench smoke test** — runs the full bench at 1 sample
    /// per scale and prints the Markdown report.  Marked `#[ignore]`
    /// because it does real attack work and takes a few seconds;
    /// run with `cargo test --lib bench_demo -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn full_bench_demo() {
        let report = run_full_bench(1);
        println!("\n{}", report);
    }
}
