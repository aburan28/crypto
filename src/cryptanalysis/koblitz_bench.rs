//! # Measurement harness for Koblitz-curve index-calculus decomposition.
//!
//! The decomposition oracle is the step that decides whether index
//! calculus on subfield curves can ever reach a curve anyone deploys.
//! [`crate::cryptanalysis::koblitz_index_calculus`] ships three of them
//! — exhaustive search, matrix-F4 over the Semaev system, and CDCL SAT
//! over the same system — and they agree on every instance tested.  What
//! was missing is *numbers*: how the instance grows with the field, how
//! the algebra behaves on it, and which of the three wins where.
//!
//! This module produces those numbers reproducibly, so the question can
//! be worked as an optimisation target rather than an opinion.  See
//! `RESEARCH_KOBLITZ_SCALING_TARGET.md` for the pre-registered
//! hypotheses these measurements are meant to settle.
//!
//! ## Two measurements, two reachable ranges
//!
//! **System structure** ([`profile_system`]) needs only the field and
//! the Frobenius-invariant subspace — no curve, no point counting, no
//! materialised factor base.  It therefore runs far past
//! [`crate::cryptanalysis::koblitz_index_calculus::MAX_N`], out to
//! `n = 63`, which is where the interesting scaling lives.  It reports
//! the unknown/equation counts, the degree, the Macaulay rank profile
//! and the **first fall degree**, using the same operational definition
//! as [`crate::cryptanalysis::ffd_harness`] so the subspace-restricted
//! numbers are directly comparable with that harness's full-field ones.
//!
//! **Oracle cost** ([`bench_instance`]) needs a real curve with a
//! prime-order subgroup, so it is capped at `MAX_N`.  It runs all three
//! oracles over the same random targets, checks they agree, and reports
//! median wall-clock per target.
//!
//! ## The cost model the target attacks
//!
//! An `m`-point decomposition chains `m − 1` copies of `S₃` over `m − 2`
//! intermediate points, and those intermediates range over the *whole*
//! field rather than the subspace:
//!
//! ```text
//!     unknowns(n, ℓ, m) = m·ℓ + (m − 2)·n.
//! ```
//!
//! The first term is what the factor base costs; the second is pure
//! bookkeeping, and it dominates.  A useful attack needs `m ≈ n/ℓ`
//! summands, so the chaining term grows like `n²/ℓ` — which is why the
//! solvable instances stop where they do, and why [`n_vars_for`] is the
//! number to drive down.
//!
//! ## What this is not
//!
//! Not a claim that any of this threatens a deployed curve.  The
//! largest instance here is `n = 63` for *structure* and `n ≤ 24` for a
//! full solve; sect163k1 is far away, and the paper's own conclusion is
//! that index calculus stays worse than rho at deployed sizes.  The
//! harness exists to find out *how far* away, with measurements instead
//! of extrapolation.

use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::binary_ecc::{BinaryPoint, F2mElement};
use crate::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, first_fall_degree, FieldStructure, MacaulayProfile, SolverEngine,
};
use crate::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, enumerate_decompose, groebner_decompose, invariant_subspace_basis,
    order_of_2_mod_n, sat_decompose, KoblitzCurve,
};

/// Unknowns in the chained `m`-point decomposition system: `m·ℓ`
/// subspace coordinates plus `(m − 2)·n` for the intermediate points.
///
/// The baseline the scaling target tries to beat.
pub fn n_vars_for(n: u32, ell: u32, m: usize) -> usize {
    m * ell as usize + m.saturating_sub(2) * n as usize
}

/// Largest `m` whose system still fits the 64-variable Boolean-monomial
/// budget, or `None` if even `m = 2` does not.
pub fn max_m_within_budget(n: u32, ell: u32, budget: usize) -> Option<usize> {
    (2..=64)
        .filter(|&m| n_vars_for(n, ell, m) <= budget)
        .next_back()
}

// ── System structure ───────────────────────────────────────────────

/// Structure of one decomposition system: size, degree, and how its
/// Macaulay matrices rank up.
#[derive(Clone, Debug)]
pub struct SystemProfile {
    /// Extension degree.
    pub n: u32,
    /// `ℓ = ord_n(2)`, the invariant subspace dimension.
    pub ell: u32,
    /// Summands.
    pub m: usize,
    /// Boolean unknowns.
    pub n_vars: usize,
    /// Boolean equations (`n` per `S₃` link).
    pub n_eqs: usize,
    /// Total degree of the system: 2 for `m = 2`, 3 once chained.
    pub degree: u32,
    /// First fall degree, by the `ffd_harness` definition: smallest
    /// `D ≥ 2` with `rank < rows` and `rank < cols`.
    pub fall_degree: Option<u32>,
    /// Per-degree Macaulay rank profile.
    pub macaulay: Vec<MacaulayProfile>,
    /// Milliseconds spent profiling.
    pub elapsed_ms: f64,
}

/// **Profile the decomposition system** for extension degree `n`, the
/// `factor_index`-th invariant subspace, and `m` summands.
///
/// The target abscissa is drawn from `x_seed`; `S₃` does not involve the
/// curve parameter `a`, and `b = 1` for every Koblitz curve, so no curve
/// is needed — only the field and the subspace.
///
/// Returns `None` when `n` has no non-trivial invariant subspace, when
/// the system exceeds the 64-variable Boolean budget, or when `n ≥ 64`.
pub fn profile_system(
    n: u32,
    factor_index: usize,
    m: usize,
    d_max: u32,
    x_seed: u64,
) -> Option<SystemProfile> {
    let (irr, basis) = invariant_subspace_basis(n, factor_index)?;
    let st = FieldStructure::new(n, &irr);
    let b = F2mElement::one(n);
    let mut rng = StdRng::seed_from_u64(x_seed);
    let x_r = F2mElement::from_biguint(&BigUint::from(rng.gen::<u64>()), n);

    let start = Instant::now();
    let sys = build_decomposition_system(&basis, &x_r, &b, m, &st)?;
    let degree = sys
        .equations
        .iter()
        .flat_map(|e| e.terms.iter())
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(0);
    let (fall_degree, macaulay) = first_fall_degree(&sys.equations, sys.n_vars, d_max);
    let elapsed_ms = start.elapsed().as_secs_f64() * 1e3;

    Some(SystemProfile {
        n,
        ell: basis.len() as u32,
        m,
        n_vars: sys.n_vars,
        n_eqs: sys.equations.len(),
        degree,
        fall_degree,
        macaulay,
        elapsed_ms,
    })
}

/// First-fall-degree statistics over several target draws.
///
/// The fall degree is a property of the *system*, and the system
/// depends on the target abscissa `x(R)`.  A single draw is therefore a
/// noisy statistic — measured on `K / F_2^9` at `m = 2`, one draw falls
/// at `D = 2` and another shows no fall to `D = 3`.  Anything built on
/// this number has to average.
#[derive(Clone, Debug)]
pub struct FfdSummary {
    /// Extension degree.
    pub n: u32,
    /// Subspace dimension.
    pub ell: u32,
    /// Summands.
    pub m: usize,
    /// Boolean unknowns.
    pub n_vars: usize,
    /// Boolean equations.
    pub n_eqs: usize,
    /// System degree.
    pub degree: u32,
    /// Draws taken.
    pub trials: usize,
    /// Smallest fall degree seen.
    pub fall_min: Option<u32>,
    /// Largest fall degree seen.
    pub fall_max: Option<u32>,
    /// Draws with no fall at or below `d_max`.
    pub no_fall: usize,
    /// Mean rank deficit `rows − rank` at `D = 2`, over the draws.
    pub mean_syzygies_d2: f64,
}

/// Profile `trials` independent target draws and summarise the fall
/// degree across them.
pub fn ffd_summary(
    n: u32,
    factor_index: usize,
    m: usize,
    d_max: u32,
    seed: u64,
    trials: usize,
) -> Option<FfdSummary> {
    let mut falls: Vec<Option<u32>> = Vec::with_capacity(trials);
    let mut syz = Vec::with_capacity(trials);
    let mut first: Option<SystemProfile> = None;
    for t in 0..trials {
        let p = profile_system(n, factor_index, m, d_max, seed.wrapping_add(t as u64))?;
        falls.push(p.fall_degree);
        if let Some(q) = p.macaulay.iter().find(|q| q.degree == 2) {
            syz.push(q.syzygies() as f64);
        }
        if first.is_none() {
            first = Some(p);
        }
    }
    let p = first?;
    let seen: Vec<u32> = falls.iter().flatten().copied().collect();
    Some(FfdSummary {
        n,
        ell: p.ell,
        m,
        n_vars: p.n_vars,
        n_eqs: p.n_eqs,
        degree: p.degree,
        trials,
        fall_min: seen.iter().copied().min(),
        fall_max: seen.iter().copied().max(),
        no_fall: falls.iter().filter(|f| f.is_none()).count(),
        mean_syzygies_d2: if syz.is_empty() {
            0.0
        } else {
            syz.iter().sum::<f64>() / syz.len() as f64
        },
    })
}

/// Summarise the fall degree for every `(n, m)` in the sweep.
pub fn sweep_ffd(
    ns: &[u32],
    ms: &[usize],
    d_max: u32,
    seed: u64,
    trials: usize,
) -> Vec<FfdSummary> {
    let mut out = Vec::new();
    for &n in ns {
        for &m in ms {
            if let Some(s) = ffd_summary(n, 0, m, d_max, seed, trials) {
                out.push(s);
            }
        }
    }
    out
}

/// Markdown table of fall-degree summaries.
pub fn format_ffd_table(rows: &[FfdSummary]) -> String {
    let mut out = String::new();
    out.push_str("| n | ℓ | m | vars | eqs | deg | FFD min | FFD max | no fall | mean syz D=2 |\n");
    out.push_str("|--:|--:|--:|-----:|----:|----:|--------:|--------:|--------:|-------------:|\n");
    for r in rows {
        let f = |v: Option<u32>| v.map(|d| d.to_string()).unwrap_or_else(|| "—".into());
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {}/{} | {:.2} |\n",
            r.n,
            r.ell,
            r.m,
            r.n_vars,
            r.n_eqs,
            r.degree,
            f(r.fall_min),
            f(r.fall_max),
            r.no_fall,
            r.trials,
            r.mean_syzygies_d2
        ));
    }
    out
}

/// Profile every `(n, m)` in the sweep, skipping what does not fit.
pub fn sweep_systems(ns: &[u32], ms: &[usize], d_max: u32, x_seed: u64) -> Vec<SystemProfile> {
    let mut out = Vec::new();
    for &n in ns {
        for &m in ms {
            if let Some(p) = profile_system(n, 0, m, d_max, x_seed) {
                out.push(p);
            }
        }
    }
    out
}

// ── Oracle cost ────────────────────────────────────────────────────

/// Timing and verdicts for one oracle over a set of targets.
#[derive(Clone, Debug)]
pub struct OracleRun {
    /// `"search"`, `"matrix-f4"` or `"sat"`.
    pub oracle: &'static str,
    /// Targets that decomposed.
    pub decomposed: usize,
    /// Targets proven not to decompose (an F4 certificate or an UNSAT).
    /// Exhaustive search never *proves* anything beyond its own sweep,
    /// so it reports its negatives here too.
    pub refuted: usize,
    /// Targets where the oracle ran out of budget — neither a
    /// decomposition nor a refutation.
    pub inconclusive: usize,
    /// Median milliseconds per target.
    pub median_ms: f64,
    /// Total milliseconds over all targets.
    pub total_ms: f64,
}

/// One benchmarked instance: all three oracles over the same targets.
#[derive(Clone, Debug)]
pub struct InstanceBench {
    /// Curve parameter.
    pub a: u8,
    /// Extension degree.
    pub n: u32,
    /// Subspace dimension.
    pub ell: u32,
    /// Summands.
    pub m: usize,
    /// Factor-base size.
    pub fb_size: usize,
    /// `π`-orbits, i.e. linear-algebra unknowns.
    pub orbits: usize,
    /// Boolean unknowns in the decomposition system.
    pub n_vars: usize,
    /// Targets tried.
    pub targets: usize,
    /// Per-oracle results.
    pub runs: Vec<OracleRun>,
    /// Targets where the three oracles did not return the same verdict.
    /// **Must be zero**; anything else is a bug, not a measurement.
    pub disagreements: usize,
}

fn median(mut xs: Vec<f64>) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = xs.len() / 2;
    if xs.len() % 2 == 0 {
        (xs[mid - 1] + xs[mid]) / 2.0
    } else {
        xs[mid]
    }
}

/// **Benchmark all three oracles** on `targets` pseudo-random points of
/// `⟨G⟩`, checking they agree.
///
/// `node_budget` bounds the F4 splitting search and `max_models` the
/// SAT model enumeration; a `None` from either with its budget spent is
/// counted as inconclusive rather than as a refutation.
///
/// Returns `None` if the curve or the factor base cannot be built.
#[allow(clippy::too_many_arguments)]
pub fn bench_instance(
    a: u8,
    n: u32,
    factor_index: usize,
    m: usize,
    targets: usize,
    seed: u64,
    node_budget: usize,
    max_models: usize,
) -> Option<InstanceBench> {
    let kc = KoblitzCurve::new(a, n)?;
    let fb = build_frobenius_factor_base(&kc, factor_index)?;
    let index_of = fb.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let g = kc.generator().clone();
    let n_vars = n_vars_for(n, fb.ell, m);

    let mut rng = StdRng::seed_from_u64(seed);
    let r_u64 = kc
        .subgroup_order
        .to_u64_digits()
        .first()
        .copied()
        .unwrap_or(2)
        .max(2);
    let points: Vec<BinaryPoint> = (0..targets)
        .map(|_| kc.mul(&g, &BigUint::from(rng.gen_range(1..r_u64))))
        .collect();

    let mut search = (0usize, 0usize, 0usize, Vec::new());
    let mut f4 = (0usize, 0usize, 0usize, Vec::new());
    let mut sat = (0usize, 0usize, 0usize, Vec::new());
    let mut disagreements = 0usize;

    for target in &points {
        let t = Instant::now();
        let by_search = enumerate_decompose(&kc, &fb, &index_of, target, m);
        search.3.push(t.elapsed().as_secs_f64() * 1e3);
        if by_search.is_some() {
            search.0 += 1;
        } else {
            search.1 += 1;
        }

        let t = Instant::now();
        let (by_f4, f4_stats) = groebner_decompose(
            &kc,
            &fb,
            &index_of,
            &st,
            target,
            m,
            SolverEngine::default(),
            node_budget,
        );
        f4.3.push(t.elapsed().as_secs_f64() * 1e3);
        match (&by_f4, f4_stats.exhausted) {
            (Some(_), _) => f4.0 += 1,
            (None, false) => f4.1 += 1,
            (None, true) => f4.2 += 1,
        }

        let t = Instant::now();
        let (by_sat, sat_stats) = sat_decompose(&kc, &fb, &index_of, &st, target, m, max_models);
        sat.3.push(t.elapsed().as_secs_f64() * 1e3);
        match (&by_sat, sat_stats.refuted) {
            (Some(_), _) => sat.0 += 1,
            (None, true) => sat.1 += 1,
            (None, false) => sat.2 += 1,
        }

        // Agreement gate: a decomposition found by one oracle must be
        // found by all, and every returned decomposition must be real.
        let verdicts = [by_search.is_some(), by_f4.is_some(), by_sat.is_some()];
        if !(verdicts.iter().all(|v| *v) || verdicts.iter().all(|v| !*v)) {
            disagreements += 1;
        }
        for found in [&by_search, &by_f4, &by_sat].into_iter().flatten() {
            let mut acc = BinaryPoint::Infinity;
            for i in found {
                acc = kc.add(&acc, &fb.points[*i]);
            }
            if &acc != target {
                disagreements += 1;
            }
        }
    }

    let pack = |name: &'static str, v: (usize, usize, usize, Vec<f64>)| OracleRun {
        oracle: name,
        decomposed: v.0,
        refuted: v.1,
        inconclusive: v.2,
        median_ms: median(v.3.clone()),
        total_ms: v.3.iter().sum(),
    };

    Some(InstanceBench {
        a,
        n,
        ell: fb.ell,
        m,
        fb_size: fb.points.len(),
        orbits: fb.orbits.len(),
        n_vars,
        targets,
        runs: vec![
            pack("search", search),
            pack("matrix-f4", f4),
            pack("sat", sat),
        ],
        disagreements,
    })
}

// ── Reporting ──────────────────────────────────────────────────────

/// Markdown table of system profiles.
pub fn format_system_table(profiles: &[SystemProfile]) -> String {
    let mut out = String::new();
    out.push_str(
        "| n | ℓ | m | vars | eqs | deg | FFD | D=2 rows×cols (rank) | D=3 rows×cols (rank) |\n",
    );
    out.push_str(
        "|--:|--:|--:|-----:|----:|----:|----:|---------------------:|---------------------:|\n",
    );
    for p in profiles {
        let cell = |d: u32| match p.macaulay.iter().find(|q| q.degree == d) {
            Some(q) => format!("{}×{} ({})", q.rows, q.cols, q.rank),
            None => "—".to_string(),
        };
        out.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            p.n,
            p.ell,
            p.m,
            p.n_vars,
            p.n_eqs,
            p.degree,
            p.fall_degree
                .map(|d| d.to_string())
                .unwrap_or_else(|| "—".into()),
            cell(2),
            cell(3),
        ));
    }
    out
}

/// Markdown table of oracle costs.
pub fn format_oracle_table(benches: &[InstanceBench]) -> String {
    let mut out = String::new();
    out.push_str(
        "| curve | n | ℓ | m | \\|F\\| | vars | oracle | found | refuted | inconc | median ms |\n",
    );
    out.push_str(
        "|:------|--:|--:|--:|------:|-----:|:-------|------:|--------:|-------:|----------:|\n",
    );
    for b in benches {
        for r in &b.runs {
            out.push_str(&format!(
                "| K_{} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.3} |\n",
                b.a,
                b.n,
                b.ell,
                b.m,
                b.fb_size,
                b.n_vars,
                r.oracle,
                r.decomposed,
                r.refuted,
                r.inconclusive,
                r.median_ms
            ));
        }
    }
    out
}

/// Machine-readable dump, for autolab runs that want to diff against a
/// previous baseline rather than read a table.
pub fn to_json(profiles: &[SystemProfile], benches: &[InstanceBench]) -> String {
    let systems: Vec<serde_json::Value> = profiles
        .iter()
        .map(|p| {
            serde_json::json!({
                "n": p.n, "ell": p.ell, "m": p.m,
                "n_vars": p.n_vars, "n_eqs": p.n_eqs, "degree": p.degree,
                "fall_degree": p.fall_degree,
                "macaulay": p.macaulay.iter().map(|q| serde_json::json!({
                    "degree": q.degree, "rows": q.rows, "cols": q.cols,
                    "rank": q.rank, "syzygies": q.syzygies(),
                })).collect::<Vec<_>>(),
                "elapsed_ms": p.elapsed_ms,
            })
        })
        .collect();
    let oracles: Vec<serde_json::Value> = benches
        .iter()
        .map(|b| {
            serde_json::json!({
                "a": b.a, "n": b.n, "ell": b.ell, "m": b.m,
                "fb_size": b.fb_size, "orbits": b.orbits, "n_vars": b.n_vars,
                "targets": b.targets, "disagreements": b.disagreements,
                "runs": b.runs.iter().map(|r| serde_json::json!({
                    "oracle": r.oracle, "decomposed": r.decomposed,
                    "refuted": r.refuted, "inconclusive": r.inconclusive,
                    "median_ms": r.median_ms, "total_ms": r.total_ms,
                })).collect::<Vec<_>>(),
            })
        })
        .collect();
    serde_json::json!({ "systems": systems, "oracles": oracles }).to_string()
}

/// The `n` values with a non-trivial Frobenius-invariant subspace of
/// dimension `≤ max_ell`, in range — the ladder every sweep walks.
pub fn subspace_ladder(n_max: u32, max_ell: u32) -> Vec<(u32, u32)> {
    (5..=n_max)
        .step_by(2)
        .filter_map(|n| order_of_2_mod_n(n).map(|l| (n, l)))
        .filter(|(n, l)| *l <= max_ell && *l < *n)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cost_model_matches_the_built_system() {
        // n_vars_for is the number the target drives down, so it had
        // better be the number the builder actually produces.
        for n in [7u32, 9, 15, 21] {
            let (_, basis) = invariant_subspace_basis(n, 0).unwrap();
            for m in [2usize, 3] {
                if let Some(p) = profile_system(n, 0, m, 2, 7) {
                    assert_eq!(
                        p.n_vars,
                        n_vars_for(n, basis.len() as u32, m),
                        "n = {n}, m = {m}"
                    );
                }
            }
        }
    }

    #[test]
    fn profiles_are_deterministic_and_structural() {
        let a = profile_system(9, 0, 2, 3, 42).unwrap();
        let b = profile_system(9, 0, 2, 3, 42).unwrap();
        assert_eq!(a.n_vars, b.n_vars);
        assert_eq!(a.fall_degree, b.fall_degree);
        assert_eq!(
            a.macaulay.iter().map(|p| p.rank).collect::<Vec<_>>(),
            b.macaulay.iter().map(|p| p.rank).collect::<Vec<_>>()
        );
        // K_a / F_2^9: ℓ = 6, so 12 unknowns, 9 equations, quadratic.
        assert_eq!((a.ell, a.n_vars, a.n_eqs, a.degree), (6, 12, 9, 2));
    }

    #[test]
    fn chaining_makes_the_system_cubic_and_much_larger() {
        let two = profile_system(9, 0, 2, 2, 1).unwrap();
        let three = profile_system(9, 0, 3, 2, 1).unwrap();
        assert_eq!(two.degree, 2);
        assert_eq!(three.degree, 3);
        // One intermediate point costs a full n unknowns.
        assert_eq!(three.n_vars - two.n_vars, 6 + 9);
    }

    #[test]
    fn the_budget_is_what_stops_the_sweep() {
        // n = 63, ℓ = 6: m = 2 fits, m = 3 needs 81 > 64 unknowns.
        assert_eq!(n_vars_for(63, 6, 2), 12);
        assert_eq!(n_vars_for(63, 6, 3), 81);
        assert_eq!(max_m_within_budget(63, 6, 64), Some(2));
        assert!(profile_system(63, 0, 3, 2, 1).is_none());
        assert!(profile_system(63, 0, 2, 2, 1).is_some());
    }

    #[test]
    fn the_ladder_is_the_n_with_small_subspaces() {
        let ladder = subspace_ladder(63, 6);
        assert!(ladder.contains(&(7, 3)));
        assert!(ladder.contains(&(15, 4)));
        assert!(ladder.contains(&(31, 5)));
        assert!(ladder.contains(&(63, 6)));
        // ord_11(2) = 10, too big for this cut.
        assert!(!ladder.iter().any(|(n, _)| *n == 11));
    }

    #[test]
    fn all_three_oracles_agree_on_a_benchmarked_instance() {
        let b = bench_instance(0, 9, 0, 2, 6, 0xB0B, 20_000, 64).unwrap();
        assert_eq!(b.disagreements, 0, "oracles disagreed: {b:?}");
        assert_eq!(b.runs.len(), 3);
        for r in &b.runs {
            assert_eq!(
                r.decomposed + r.refuted + r.inconclusive,
                b.targets,
                "{} lost a target",
                r.oracle
            );
        }
        // Every target of K_0/F_2^9 decomposes into 2 factor-base points.
        assert!(b.runs.iter().all(|r| r.decomposed == b.targets));
    }

    #[test]
    fn fall_degree_varies_with_the_target_draw() {
        // The reason the harness averages: one draw is not the system's
        // fall degree, it is one system's fall degree.
        let s = ffd_summary(9, 0, 2, 3, 0x5EED, 16).unwrap();
        assert_eq!(s.trials, 16);
        assert!(
            s.no_fall > 0 || s.fall_min != s.fall_max,
            "expected spread across draws, got {s:?}"
        );
        // …and the summary still describes the same system.
        assert_eq!((s.ell, s.n_vars, s.n_eqs, s.degree), (6, 12, 9, 2));
    }

    #[test]
    fn json_round_trips_as_valid_json() {
        let profiles = vec![profile_system(9, 0, 2, 3, 1).unwrap()];
        let benches = vec![bench_instance(0, 9, 0, 2, 2, 1, 20_000, 64).unwrap()];
        let text = to_json(&profiles, &benches);
        let parsed: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(parsed["systems"][0]["n_vars"], 12);
        assert_eq!(parsed["oracles"][0]["disagreements"], 0);
    }
}
