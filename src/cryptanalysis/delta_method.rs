//! # The Shantz–Teske "delta method" for Semaev decomposition systems.
//!
//! Shantz & Teske (*Solving the elliptic curve discrete logarithm problem
//! using Semaev polynomials, Weil descent and Gröbner basis methods — an
//! experimental study*, 2013/2015) observe that with a subspace factor
//! base of dimension `ℓ` and `m` summands the Weil-descended system has
//! `m·ℓ` unknowns and `n` equations, and that choosing
//!
//! ```text
//!     m·ℓ = n − δ,   δ > 0
//! ```
//!
//! makes it overdetermined: each system is solved faster (lower solving
//! degree, fewer spurious solutions) at the price of a `2^{−δ}` lower
//! decomposition probability per target.  The crate's ladder already
//! parameterises `ℓ` but nothing reports the **net** effect, which is
//! what decides collection:
//!
//! ```text
//!     cost per relation = (1 / P[target decomposes]) × (wall per system)
//! ```
//!
//! [`delta_dimensions`] enumerates the admissible `ℓ` for a `δ` window,
//! [`predicted_decomposition_probability`] gives the subspace coverage
//! model the measurement is compared against, and [`run_delta_cell`]
//! measures one `(n, m, ℓ)` cell on seeded public targets through the
//! existing [`groebner_decompose`] oracle.  `examples/delta_method_sweep.rs`
//! sweeps and reports.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};

use crate::binary_ecc::BinaryPoint;

use super::koblitz_groebner::{f4_word_ops_thread, FieldStructure, SolverEngine, MAX_VARS};
use super::koblitz_index_calculus::{
    build_standard_subspace_factor_base, groebner_decompose, point_key, FrobeniusFactorBase,
    KoblitzCurve,
};

/// Largest subspace dimension [`build_standard_subspace_factor_base`]
/// accepts.
pub const MAX_STANDARD_DIMENSION: u32 = 20;

/// One admissible dimension for a `δ` window.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DeltaChoice {
    /// Subspace dimension.
    pub ell: u32,
    /// `δ = n − m·ℓ`: positive is overdetermined, negative underdetermined.
    pub delta: i32,
    /// Boolean unknowns of the chained system, `m·ℓ + (m − 2)·n`.
    pub n_vars: usize,
}

/// Every `ℓ` with `n − m·ℓ` inside `[delta_min, delta_max]`, restricted to
/// what the system builder can hold (`1 ≤ ℓ ≤ min(20, n − 1)`, at most
/// [`MAX_VARS`] unknowns), in increasing `ℓ` (decreasing `δ`).  A `δ` whose
/// `(n − δ)` is not a multiple of `m` has no exact `ℓ`; the window is what
/// makes "clamped to valid `ℓ`" well defined.
pub fn delta_dimensions(n: u32, m: usize, delta_min: i32, delta_max: i32) -> Vec<DeltaChoice> {
    let mut out = Vec::new();
    if m < 2 || n < 3 || delta_min > delta_max {
        return out;
    }
    let max_ell = MAX_STANDARD_DIMENSION.min(n - 1);
    for ell in 1..=max_ell {
        let delta = n as i64 - (m as i64) * (ell as i64);
        if delta < delta_min as i64 || delta > delta_max as i64 {
            continue;
        }
        let n_vars = m * ell as usize + (m - 2) * n as usize;
        if n_vars > MAX_VARS {
            continue;
        }
        out.push(DeltaChoice {
            ell,
            delta: delta as i32,
            n_vars,
        });
    }
    out
}

/// `2^{−δ}`: the relative probability factor Shantz–Teske predict.
pub fn delta_probability_factor(delta: i32) -> f64 {
    (2f64).powi(-delta)
}

/// Subspace coverage model: a random target of a group of order
/// `group_order` is an `m`-sum of a base of `points` points with
/// probability `1 − exp(−|F|^m / (m!·#E))` — the `|F|^m/m!` unordered
/// sums land uniformly on the `#E` group elements, and only the
/// fraction `1/h` that lands in the subgroup can hit a subgroup target,
/// which is why `#E` and not `r` divides.
pub fn predicted_decomposition_probability(points: usize, m: usize, group_order: f64) -> f64 {
    let factorial: f64 = (1..=m).map(|k| k as f64).product();
    let expected = (points as f64).powi(m as i32) / (factorial * group_order);
    1.0 - (-expected).exp()
}

/// One target's measurement.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DeltaTargetRecord {
    pub scalar: u64,
    pub decomposed: bool,
    pub wall_ms: f64,
    pub word_xors: u64,
    pub max_degree_built: u32,
    pub splits: usize,
    pub reductions: usize,
    pub exhausted: bool,
}

/// One `(n, m, ℓ)` cell.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DeltaCellResult {
    pub n: u32,
    pub a: u8,
    pub m: usize,
    pub ell: u32,
    pub delta: i32,
    pub n_vars: usize,
    pub base_points: usize,
    pub group_order: String,
    pub subgroup_order: String,
    pub targets_requested: usize,
    pub targets_measured: usize,
    pub decomposed: usize,
    /// `decomposed / targets_measured`.
    pub measured_fraction: f64,
    /// Coverage model with the actual `|F|`.
    pub predicted_fraction: f64,
    /// `2^{−δ}`.
    pub predicted_delta_factor: f64,
    pub mean_wall_ms: f64,
    pub median_wall_ms: f64,
    pub mean_word_xors: f64,
    /// Highest Macaulay degree built on any target, and whether every
    /// target was decided without splitting.
    pub max_degree_built: u32,
    pub split_free: bool,
    pub mean_splits: f64,
    /// `mean_wall_ms / measured_fraction`; infinite when nothing
    /// decomposed.
    pub net_ms_per_relation: f64,
    pub net_word_xors_per_relation: f64,
    /// The cell stopped at its deadline before every target ran.
    pub truncated: bool,
    pub records: Vec<DeltaTargetRecord>,
}

/// Seeded public targets `[k]G` with their scalars (splitmix64 stream).
pub fn seeded_targets(kc: &KoblitzCurve, count: usize, seed: u64) -> Vec<(u64, BinaryPoint)> {
    let r = kc.subgroup_order.to_u64().unwrap_or(u64::MAX);
    let mut state = seed ^ 0xde17_a000_0000_5e7a;
    let g = kc.generator();
    let mut out = Vec::with_capacity(count);
    while out.len() < count {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        let k = 1 + z % (r - 1);
        out.push((k, kc.mul(g, &BigUint::from(k))));
    }
    out
}

/// The standard subspace base of dimension `ell` and its index map.
pub fn standard_base(
    kc: &KoblitzCurve,
    ell: u32,
) -> Result<(FrobeniusFactorBase, HashMap<(BigUint, BigUint), usize>), String> {
    let base = build_standard_subspace_factor_base(kc, ell)?;
    let index_of = base
        .points
        .iter()
        .enumerate()
        .map(|(i, p)| (point_key(p), i))
        .collect();
    Ok((base, index_of))
}

/// Measure one cell: every target through [`groebner_decompose`] on the
/// standard subspace of dimension `ell`, until `deadline` if given.
#[allow(clippy::too_many_arguments)]
pub fn run_delta_cell(
    kc: &KoblitzCurve,
    m: usize,
    ell: u32,
    targets: &[(u64, BinaryPoint)],
    engine: SolverEngine,
    node_budget: usize,
    deadline: Option<Instant>,
    mut on_record: impl FnMut(&DeltaTargetRecord),
) -> Result<DeltaCellResult, String> {
    let (base, index_of) = standard_base(kc, ell)?;
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let n_vars = m * ell as usize + (m - 2) * kc.n as usize;
    if n_vars > MAX_VARS {
        return Err(format!(
            "{n_vars} unknowns exceed the {MAX_VARS}-variable engine"
        ));
    }
    let mut records = Vec::with_capacity(targets.len());
    let mut truncated = false;
    for (scalar, target) in targets {
        if deadline.is_some_and(|d| Instant::now() >= d) {
            truncated = true;
            break;
        }
        let before = f4_word_ops_thread();
        let started = Instant::now();
        let (found, stats) =
            groebner_decompose(kc, &base, &index_of, &st, target, m, engine, node_budget);
        let wall_ms = started.elapsed().as_secs_f64() * 1e3;
        let record = DeltaTargetRecord {
            scalar: *scalar,
            decomposed: found.is_some(),
            wall_ms,
            word_xors: f4_word_ops_thread().saturating_sub(before),
            max_degree_built: stats.max_degree_built,
            splits: stats.splits,
            reductions: stats.reductions,
            exhausted: stats.exhausted,
        };
        on_record(&record);
        records.push(record);
    }
    Ok(summarise_cell(
        kc,
        m,
        ell,
        &base,
        targets.len(),
        records,
        truncated,
    ))
}

/// Fold the records of a cell into its summary.
pub fn summarise_cell(
    kc: &KoblitzCurve,
    m: usize,
    ell: u32,
    base: &FrobeniusFactorBase,
    targets_requested: usize,
    records: Vec<DeltaTargetRecord>,
    truncated: bool,
) -> DeltaCellResult {
    let measured = records.len();
    let decomposed = records.iter().filter(|r| r.decomposed).count();
    let fraction = if measured == 0 {
        0.0
    } else {
        decomposed as f64 / measured as f64
    };
    let mean_wall = if measured == 0 {
        0.0
    } else {
        records.iter().map(|r| r.wall_ms).sum::<f64>() / measured as f64
    };
    let mut walls: Vec<f64> = records.iter().map(|r| r.wall_ms).collect();
    walls.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_wall = if walls.is_empty() {
        0.0
    } else {
        walls[walls.len() / 2]
    };
    let mean_xors = if measured == 0 {
        0.0
    } else {
        records.iter().map(|r| r.word_xors as f64).sum::<f64>() / measured as f64
    };
    let mean_splits = if measured == 0 {
        0.0
    } else {
        records.iter().map(|r| r.splits as f64).sum::<f64>() / measured as f64
    };
    let delta = kc.n as i32 - (m as i32) * ell as i32;
    let group_order = kc.group_order.to_f64().unwrap_or(f64::INFINITY);
    DeltaCellResult {
        n: kc.n,
        a: kc.a,
        m,
        ell,
        delta,
        n_vars: m * ell as usize + (m - 2) * kc.n as usize,
        base_points: base.points.len(),
        group_order: kc.group_order.to_string(),
        subgroup_order: kc.subgroup_order.to_string(),
        targets_requested,
        targets_measured: measured,
        decomposed,
        measured_fraction: fraction,
        predicted_fraction: predicted_decomposition_probability(base.points.len(), m, group_order),
        predicted_delta_factor: delta_probability_factor(delta),
        mean_wall_ms: mean_wall,
        median_wall_ms: median_wall,
        mean_word_xors: mean_xors,
        max_degree_built: records
            .iter()
            .map(|r| r.max_degree_built)
            .max()
            .unwrap_or(0),
        split_free: records.iter().all(|r| r.splits == 0),
        mean_splits,
        net_ms_per_relation: if fraction > 0.0 {
            mean_wall / fraction
        } else {
            f64::INFINITY
        },
        net_word_xors_per_relation: if fraction > 0.0 {
            mean_xors / fraction
        } else {
            f64::INFINITY
        },
        truncated,
        records,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_index_calculus::PairSumTable;

    #[test]
    fn dimension_selection_follows_m_ell_equals_n_minus_delta() {
        let n17_m2 = delta_dimensions(17, 2, -2, 4);
        assert_eq!(
            n17_m2.iter().map(|c| (c.ell, c.delta)).collect::<Vec<_>>(),
            vec![(7, 3), (8, 1), (9, -1)]
        );
        let n17_m3 = delta_dimensions(17, 3, -2, 4);
        assert_eq!(
            n17_m3
                .iter()
                .map(|c| (c.ell, c.delta, c.n_vars))
                .collect::<Vec<_>>(),
            vec![(5, 2, 32), (6, -1, 35)]
        );
        let n31_m3 = delta_dimensions(31, 3, -2, 4);
        assert_eq!(
            n31_m3
                .iter()
                .map(|c| (c.ell, c.delta, c.n_vars))
                .collect::<Vec<_>>(),
            vec![(9, 4, 58), (10, 1, 61), (11, -2, 64)]
        );
        // The variable cap clamps: at n = 41, m = 3 every ℓ with δ ≤ 4
        // needs more than 64 unknowns.
        assert!(delta_dimensions(41, 3, -2, 4).is_empty());
        // The dimension cap clamps too: n = 61, m = 2 asks for ℓ ≥ 29.
        assert!(delta_dimensions(61, 2, -2, 4).is_empty());
        for c in delta_dimensions(23, 2, -2, 4) {
            assert_eq!(c.delta, 23 - 2 * c.ell as i32);
            assert!(c.ell <= MAX_STANDARD_DIMENSION);
        }
        assert!(delta_dimensions(17, 1, -2, 4).is_empty());
        assert!(delta_dimensions(17, 2, 4, -2).is_empty());
        assert_eq!(delta_probability_factor(2), 0.25);
        assert_eq!(delta_probability_factor(-1), 2.0);
    }

    /// At a small cell the solver's yield fraction equals the exact
    /// pair-table fraction and sits within the coverage model's band.
    #[test]
    fn yield_fraction_matches_pair_table_and_coverage_model_at_a_small_cell() {
        // K_1: K_0 at n = 11 has #E = 4·23² and no usable subgroup.
        let kc = KoblitzCurve::new(1, 11).expect("K_1 / F_2^11");
        let (m, ell) = (2usize, 5u32);
        let choice = delta_dimensions(11, m, -2, 4)
            .into_iter()
            .find(|c| c.ell == ell)
            .expect("ℓ = 5 is admissible at n = 11, m = 2");
        assert_eq!(choice.delta, 1);
        let targets = seeded_targets(&kc, 48, 11);
        let cell = run_delta_cell(
            &kc,
            m,
            ell,
            &targets,
            SolverEngine::MatrixF4 { max_degree: 3 },
            4096,
            None,
            |_| {},
        )
        .expect("cell runs");
        assert_eq!(cell.targets_measured, 48);
        assert!(!cell.truncated);
        // Exact ground truth from the pair table on the same base.
        let (base, _) = standard_base(&kc, ell).expect("base");
        let table = PairSumTable::build_full(&kc, &base).expect("pair table");
        let exact: Vec<bool> = targets
            .iter()
            .map(|(_, t)| table.decompose(&kc, &base, t, m).is_some())
            .collect();
        let exact_count = exact.iter().filter(|&&d| d).count();
        assert_eq!(
            cell.decomposed, exact_count,
            "solver and pair table disagree"
        );
        for (record, truth) in cell.records.iter().zip(&exact) {
            assert_eq!(record.decomposed, *truth, "scalar {}", record.scalar);
        }
        // The model: |F|²/(2·#E) expected witnesses per target.
        assert!(cell.predicted_fraction > 0.05 && cell.predicted_fraction < 0.95);
        let ratio = cell.measured_fraction / cell.predicted_fraction;
        assert!(
            (0.4..=2.5).contains(&ratio),
            "measured {} vs predicted {} (ratio {ratio})",
            cell.measured_fraction,
            cell.predicted_fraction
        );
        assert!(cell.net_ms_per_relation.is_finite());
    }
}
