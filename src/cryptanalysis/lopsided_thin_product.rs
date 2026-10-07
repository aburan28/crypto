//! # Lopsided thin products and All-Edges Sparse Triangle: IC integration record.
//!
//! This module integrates the thin-matrix-product technique of Alman and
//! Vassilevska Williams (`arXiv:2610.06783v1`, 5 Oct 2026) into this
//! repository's index-calculus (IC) pipeline as a **stage diagnostic with no
//! end-to-end speedup claim**.
//!
//! ## What the paper shows (derived, not measured here)
//!
//! Let `X` be an `N x D` integer matrix and `Y` a `D x N` integer matrix with
//! `D <= N^{1/18}`, and let `W` be any set of at most `N^2 / sqrt(D)`
//! positions. The paper's Theorem 1 (Corollary 26 / Theorem 25) computes the
//! entries `(X Y)[i,j]` for `(i,j)` in `W` deterministically in
//! `O(N^2 / D^0.063)` operations on `O(log N)`-bit integers. That is
//! polynomially less than writing down `X Y` or than computing `|W|`
//! inner products one by one (`N^2 sqrt(D)`). More generally, for every
//! `epsilon < 0.1204` and every `kappa > 0` there is a `gamma > 0` with cost
//! `O(N^2 / D^gamma)` whenever `D <= N^epsilon` and `|W| <= N^2 / D^kappa`.
//!
//! In graph language this is the counting version of **Lopsided All-Edges
//! Sparse Triangle**: in a tripartite graph with two parts of `n` vertices
//! and a middle part of `n^epsilon` vertices, with `X` and `Y` the two
//! biadjacency matrices, `(X Y)[a,b]` is the number of middle vertices
//! adjacent to both `a` and `b`. For `epsilon < 1/18` the paper counts `|W|`
//! prescribed outer pairs in
//! `O(|W| n^{0.437 epsilon} + n^{2 - 0.063 epsilon})` time, which is
//! `O(n^{2 - 0.063 epsilon})` for `|W| <= n^{2 - epsilon/2}`. The data
//! structure version (Theorem 3 / Theorem 24) preprocesses `(X, Y)`
//! deterministically in `O(N^2 / D^0.063)` time and answers any single entry
//! in `O(D^0.437)` time.
//!
//! The technique modifies a variant of Coppersmith's rectangular matrix
//! multiplication algorithm, built from Schoenhage's ten-multiplication
//! identity, to visit only the recursion-tree leaves that the wanted entries
//! in `W` need, sharing encoding work across blocks. By known reductions
//! (Patrascu; Vassilevska Williams-Xu; Chan-Xu; Chan-Vassilevska
//! Williams-Xu), Exact Triangle, 3SUM, and APSP reduce to this lopsided
//! problem, which is how the paper refutes the 3SUM and APSP hypotheses.
//!
//! ## What this module does
//!
//! - [`ThinProductSpec`]: the exact problem statement (`N`, `D`, wanted set
//!   `W`) with validation and regime predicates (`D^18 <= N`,
//!   `|W| <= N^2 / sqrt(D)`, `epsilon`/`kappa` bookkeeping).
//! - [`wanted_product_i64`]: a checked reference implementation that computes
//!   exactly the wanted entries by direct inner products. This is the
//!   correctness baseline the paper's faster algorithm must match; it is
//!   **not** the paper's algorithm.
//! - [`paper_ops_bound`] / [`paper_query_bound`] / [`graph_time_bound`]:
//!   the paper's quoted asymptotic bounds as evaluable `f64` cost-model
//!   hooks, clearly labelled as quoted bounds rather than measurements.
//! - [`count_triangles_through_pairs`]: the lopsided-graph reading of the
//!   same computation over Boolean biadjacency matrices.
//! - [`IcInsertionPoint`] / [`insertion_assessment`]: the four IC pipeline
//!   points where a thin-product primitive could in principle be tried,
//!   each assessed as a **proposal** with the frozen experiment that would be
//!   needed to promote it. None is wired into the end-to-end `S` accounting.
//!
//! ## What this module does not do
//!
//! It does not implement the paper's pruned-recursion algorithm, does not
//! change any IC candidate's `S = total operations / sqrt(r)` cost, does not
//! re-match any rho reference, and does not update the scoreboard,
//! leaderboard, or lab browser. A row without a verified answer is not a
//! result (`AGENTS.md` section 2); this module produces no such row, so every
//! panel that cites IC/rho ratios remains current without modification. The
//! accompanying study (`research/lopsided_thin_product_20261006/`) records
//! the graphs checked and why they are unchanged.
//!
//! ## Boundary, table, ratio
//!
//! The applicable boundaries are unchanged: the generic-group floor and the
//! matched rho reference measured in the same unit on the same instances
//! (`AGENTS.md` sections 1-2). A future thin-product IC variant would earn a
//! table row only as a full-pipeline `S` measurement with every phase priced
//! (setup, queries, PDP, verification, relation matrix, linear algebra,
//! descent, recovery check) against those boundaries, classified per section
//! 3 as advance, engineering, relabelling, or accounting. Reducing a stage
//! count (wanted-entry inner products, triangle counts, query time) without
//! moving total `S` is the relabelling class, not a speedup.

use std::collections::BTreeSet;

/// Exact statement of one thin-product task: compute `(X Y)[i,j]` for
/// `(i,j)` in `wanted`, where `X` is `n_outer x middle_dim` and `Y` is
/// `middle_dim x n_outer`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ThinProductSpec {
    /// `N`: number of rows of `X` / columns of `Y` (outer dimension).
    pub n_outer: usize,
    /// `D`: inner (middle) dimension.
    pub middle_dim: usize,
    /// Wanted output positions, each `(row, col)` with both `< n_outer`.
    pub wanted: Vec<(usize, usize)>,
}

/// Validation failures for [`ThinProductSpec`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpecError {
    /// `n_outer == 0`.
    EmptyOuter,
    /// `middle_dim == 0`.
    EmptyMiddle,
    /// A wanted position lies outside the `N x N` output.
    WantedOutOfRange {
        row: usize,
        col: usize,
        n_outer: usize,
    },
    /// The wanted set contains a duplicate position.
    DuplicateWanted { row: usize, col: usize },
}

impl std::fmt::Display for SpecError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpecError::EmptyOuter => write!(f, "n_outer must be nonzero"),
            SpecError::EmptyMiddle => write!(f, "middle_dim must be nonzero"),
            SpecError::WantedOutOfRange { row, col, n_outer } => write!(
                f,
                "wanted position ({row},{col}) outside {n_outer}x{n_outer} output"
            ),
            SpecError::DuplicateWanted { row, col } => {
                write!(f, "duplicate wanted position ({row},{col})")
            }
        }
    }
}

impl std::error::Error for SpecError {}

impl ThinProductSpec {
    /// Validate dimensions and the wanted set (range + uniqueness).
    pub fn validate(&self) -> Result<(), SpecError> {
        if self.n_outer == 0 {
            return Err(SpecError::EmptyOuter);
        }
        if self.middle_dim == 0 {
            return Err(SpecError::EmptyMiddle);
        }
        let mut seen = BTreeSet::new();
        for &(row, col) in &self.wanted {
            if row >= self.n_outer || col >= self.n_outer {
                return Err(SpecError::WantedOutOfRange {
                    row,
                    col,
                    n_outer: self.n_outer,
                });
            }
            if !seen.insert((row, col)) {
                return Err(SpecError::DuplicateWanted { row, col });
            }
        }
        Ok(())
    }

    /// Number of wanted entries `|W|`.
    pub fn wanted_count(&self) -> usize {
        self.wanted.len()
    }

    /// `epsilon = ln D / ln N`, the lopsided exponent (`D = N^epsilon`).
    /// Returns `None` for degenerate `N < 2` or `D < 1`.
    pub fn epsilon(&self) -> Option<f64> {
        if self.n_outer < 2 || self.middle_dim < 1 {
            return None;
        }
        Some((self.middle_dim as f64).ln() / (self.n_outer as f64).ln())
    }

    /// `kappa` defined by `|W| = N^2 / D^kappa`, i.e.
    /// `kappa = ln(N^2/|W|) / ln D`. Returns `None` when undefined
    /// (empty `W`, `D < 2`).
    pub fn kappa(&self) -> Option<f64> {
        if self.wanted.is_empty() || self.middle_dim < 2 {
            return None;
        }
        let n2 = (self.n_outer as f64).powi(2);
        let w = self.wanted.len() as f64;
        if w <= 0.0 || w > n2 {
            return None;
        }
        Some((n2 / w).ln() / (self.middle_dim as f64).ln())
    }

    /// Paper's concrete regime: `N >= D^18` (equivalently
    /// `epsilon <= 1/18`), the setting of Theorem 5 / Theorem 1.
    pub fn in_concrete_regime(&self) -> bool {
        let d = self.middle_dim as f64;
        let n = self.n_outer as f64;
        if n < 1.0 || d < 1.0 {
            return false;
        }
        n >= d.powf(18.0)
    }

    /// Paper's general-regime check: `D <= N^epsilon_max`.
    pub fn in_epsilon_regime(&self, epsilon_max: f64) -> bool {
        match self.epsilon() {
            Some(e) => e <= epsilon_max,
            None => false,
        }
    }

    /// Wanted-set size check from the paper: `|W| <= N^2 / sqrt(D)`,
    /// i.e. `kappa >= 1/2`. An empty wanted set trivially satisfies it.
    pub fn wanted_within_paper_bound(&self) -> bool {
        if self.wanted.is_empty() {
            return true;
        }
        let d = self.middle_dim as f64;
        let n = self.n_outer as f64;
        (self.wanted.len() as f64) <= n * n / d.sqrt()
    }

    /// Both gates of the paper's concrete theorem: `N >= D^18` and
    /// `|W| <= N^2 / sqrt(D)`.
    pub fn satisfies_concrete_theorem_gates(&self) -> bool {
        self.in_concrete_regime() && self.wanted_within_paper_bound()
    }
}

/// Compute exactly the wanted entries of `X Y` by direct inner products.
///
/// `x` is `N x D` and `y` is `D x N` in row-major order. Accumulates in
/// `i128` and rejects overflow of the `i64` output. This reference
/// implementation costs `|W| * D` scalar multiplications; it establishes
/// what the paper's faster algorithm must return, not how fast it runs.
pub fn wanted_product_i64(
    spec: &ThinProductSpec,
    x: &[i64],
    y: &[i64],
) -> Result<Vec<i64>, String> {
    spec.validate().map_err(|e| e.to_string())?;
    let n = spec.n_outer;
    let d = spec.middle_dim;
    if x.len() != n * d {
        return Err(format!("x has {} entries, need N*D = {}", x.len(), n * d));
    }
    if y.len() != d * n {
        return Err(format!("y has {} entries, need D*N = {}", y.len(), d * n));
    }
    let mut out = Vec::with_capacity(spec.wanted.len());
    for &(i, j) in &spec.wanted {
        let mut acc: i128 = 0;
        for k in 0..d {
            acc += (x[i * d + k] as i128) * (y[k * n + j] as i128);
        }
        if acc < i64::MIN as i128 || acc > i64::MAX as i128 {
            return Err(format!("inner product ({i},{j}) overflows i64"));
        }
        out.push(acc as i64);
    }
    Ok(out)
}

/// Paper's concrete preprocessing/offline bound as a `f64` hook:
/// `N^2 / D^0.063` (Theorem 1). Quoted bound, not a measurement.
pub fn paper_ops_bound(n_outer: usize, middle_dim: usize) -> Option<f64> {
    if n_outer == 0 || middle_dim == 0 {
        return None;
    }
    let n = n_outer as f64;
    let d = middle_dim as f64;
    Some(n * n / d.powf(0.063))
}

/// Paper's data-structure query bound as a `f64` hook: `D^0.437`
/// (Theorem 3). Quoted bound, not a measurement.
pub fn paper_query_bound(middle_dim: usize) -> Option<f64> {
    if middle_dim == 0 {
        return None;
    }
    Some((middle_dim as f64).powf(0.437))
}

/// Lopsided-graph form of the bound for `epsilon < 1/18`:
/// `|W| n^{0.437 epsilon} + n^{2 - 0.063 epsilon}` (Section 1.1).
/// Quoted bound, not a measurement.
pub fn graph_time_bound(n: usize, epsilon: f64, wanted: usize) -> Option<f64> {
    if n == 0 || epsilon <= 0.0 {
        return None;
    }
    let nf = n as f64;
    let wf = wanted as f64;
    Some(wf * nf.powf(0.437 * epsilon) + nf.powf(2.0 - 0.063 * epsilon))
}

/// Naive per-entry baseline the paper beats in its regime: `|W| * D`
/// scalar multiplications (one inner product per wanted entry).
pub fn naive_inner_product_cost(spec: &ThinProductSpec) -> Option<u64> {
    (spec.wanted_count() as u64).checked_mul(spec.middle_dim as u64)
}

/// Counting version of Lopsided All-Edges Sparse Triangle over Boolean
/// biadjacency matrices: `x` (`N x D`) and `y` (`D x N`) hold 0/1 entries;
/// each wanted outer pair `(a,b)` counts middle vertices adjacent to both.
/// Equals the integer wanted product on 0/1 inputs.
pub fn count_triangles_through_pairs(
    spec: &ThinProductSpec,
    x01: &[u8],
    y01: &[u8],
) -> Result<Vec<u64>, String> {
    spec.validate().map_err(|e| e.to_string())?;
    let n = spec.n_outer;
    let d = spec.middle_dim;
    if x01.len() != n * d {
        return Err(format!(
            "x01 has {} entries, need N*D = {}",
            x01.len(),
            n * d
        ));
    }
    if y01.len() != d * n {
        return Err(format!(
            "y01 has {} entries, need D*N = {}",
            y01.len(),
            d * n
        ));
    }
    for &b in x01.iter().chain(y01.iter()) {
        if b > 1 {
            return Err("biadjacency matrices must be 0/1".to_string());
        }
    }
    let mut out = Vec::with_capacity(spec.wanted.len());
    for &(a, b) in &spec.wanted {
        let mut count: u64 = 0;
        for k in 0..d {
            if x01[a * d + k] == 1 && y01[k * n + b] == 1 {
                count += 1;
            }
        }
        out.push(count);
    }
    Ok(out)
}

/// Candidate IC pipeline points where a thin-product primitive could be
/// tried. All variants are proposals; none is implemented or measured here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IcInsertionPoint {
    /// Batch-verify many candidate PDP decompositions at once.
    RelationVerificationBatch,
    /// Answer many factor-base membership/intersection queries, offline
    /// (`W` known in advance) or online (data-structure version).
    FactorBaseIntersectionQuery,
    /// Accelerate one block step of sparse linear algebra (Wiedemann/Lanczos).
    SparseLinearAlgebraBlock,
    /// Batch the recursive decompositions inside target descent.
    TargetDescentBatch,
}

/// Assessment of one insertion point: why it does not inherit the paper's
/// bound and what frozen experiment would be needed to promote it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InsertionAssessment {
    /// Which pipeline point was assessed.
    pub point: IcInsertionPoint,
    /// Whether any code path in this repository currently uses a thin
    /// product at this point.
    pub implemented: bool,
    /// Short reason the paper's bound does not transfer as stated.
    pub blocking_difference: &'static str,
    /// Frozen experiment that would be needed for an end-to-end claim.
    pub required_experiment: &'static str,
}

/// Assess one IC insertion point. Every arm returns `implemented: false`:
/// recording the proposal is the deliverable; wiring it into `S` accounting
/// requires a matched baseline/candidate full-DLP comparison that does not
/// exist yet.
pub fn insertion_assessment(point: IcInsertionPoint) -> InsertionAssessment {
    match point {
        IcInsertionPoint::RelationVerificationBatch => InsertionAssessment {
            point,
            implemented: false,
            blocking_difference: "PDP systems are polynomial (Semaev/Weil descent), not inner products; batching them as dense thin products needs a new reduction, not a substitution",
            required_experiment: "freeze curve/subgroup/base/targets/seeds; run matched baseline/candidate full-DLP suites with identical inputs and price every phase in S; stage-only PDP counts stay a diagnostic",
        },
        IcInsertionPoint::FactorBaseIntersectionQuery => InsertionAssessment {
            point,
            implemented: false,
            blocking_difference: "closest fit (offline W vs online queries mirrors relation queries), but IC yield depends on algebraic decomposition probability, not just set intersection; swapping the oracle changes the method and its floor ratio",
            required_experiment: "freeze one-target workloads; pair candidate IC against matched rho on the same public targets under the same envelope; report rho_online/ic_online with all five online phases named",
        },
        IcInsertionPoint::SparseLinearAlgebraBlock => InsertionAssessment {
            point,
            implemented: false,
            blocking_difference: "block Wiedemann needs large sparse matrix-vector products over Z/rZ, not small-integer dense thin products; shapes and rings differ",
            required_experiment: "price relation-matrix build plus final LA inside total S on the same relation sets; a kernel-only microbenchmark cannot promote a full-pipeline claim",
        },
        IcInsertionPoint::TargetDescentBatch => InsertionAssessment {
            point,
            implemented: false,
            blocking_difference: "descent is recursive per-target decomposition work charged once per target; batching across targets changes the workload from one-target to multi-target and cannot replace the primary comparison",
            required_experiment: "complete the single-target measurement first, then pose any multi-target question separately with shared setup labelled and k named",
        },
    }
}

/// All four insertion points in a fixed order, for reports and tables.
pub fn all_insertion_points() -> [IcInsertionPoint; 4] {
    [
        IcInsertionPoint::RelationVerificationBatch,
        IcInsertionPoint::FactorBaseIntersectionQuery,
        IcInsertionPoint::SparseLinearAlgebraBlock,
        IcInsertionPoint::TargetDescentBatch,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny_spec() -> ThinProductSpec {
        ThinProductSpec {
            n_outer: 4,
            middle_dim: 2,
            wanted: vec![(0, 0), (0, 3), (2, 1)],
        }
    }

    #[test]
    fn validate_accepts_wellformed_spec() {
        assert!(tiny_spec().validate().is_ok());
        assert_eq!(tiny_spec().wanted_count(), 3);
    }

    #[test]
    fn validate_rejects_bad_specs() {
        assert_eq!(
            ThinProductSpec {
                n_outer: 0,
                middle_dim: 2,
                wanted: vec![]
            }
            .validate(),
            Err(SpecError::EmptyOuter)
        );
        assert_eq!(
            ThinProductSpec {
                n_outer: 3,
                middle_dim: 0,
                wanted: vec![]
            }
            .validate(),
            Err(SpecError::EmptyMiddle)
        );
        assert!(matches!(
            ThinProductSpec {
                n_outer: 2,
                middle_dim: 2,
                wanted: vec![(2, 0)]
            }
            .validate(),
            Err(SpecError::WantedOutOfRange { .. })
        ));
        assert!(matches!(
            ThinProductSpec {
                n_outer: 2,
                middle_dim: 2,
                wanted: vec![(0, 0), (0, 0)]
            }
            .validate(),
            Err(SpecError::DuplicateWanted { .. })
        ));
    }

    #[test]
    fn wanted_product_matches_naive_full_product() {
        // X = [[1,2],[3,4],[5,6],[7,8]], Y = [[1,0,1,0],[0,1,1,1]].
        let spec = tiny_spec();
        let x: Vec<i64> = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let y: Vec<i64> = vec![1, 0, 1, 0, 0, 1, 1, 1];
        let got = wanted_product_i64(&spec, &x, &y).unwrap();
        // (0,0)=1, (0,3)=2, (2,1)=6.
        assert_eq!(got, vec![1, 2, 6]);
    }

    #[test]
    fn wanted_product_rejects_shape_mismatch_and_overflow() {
        let spec = tiny_spec();
        assert!(wanted_product_i64(&spec, &[1, 2], &[1; 8]).is_err());
        assert!(wanted_product_i64(&spec, &[1; 8], &[1; 2]).is_err());
        let big_x = vec![i64::MAX; 8];
        let big_y = vec![2_i64; 8];
        assert!(wanted_product_i64(&spec, &big_x, &big_y).is_err());
    }

    #[test]
    fn regime_gates_behave() {
        // D=4, N=4^18 satisfies N >= D^18 exactly; use smaller analogue
        // D=2, N=2^18 for the predicate shape.
        let spec = ThinProductSpec {
            n_outer: 262144, // 2^18
            middle_dim: 2,
            wanted: vec![(0, 0)],
        };
        assert!(spec.in_concrete_regime());
        assert!(spec.wanted_within_paper_bound());
        assert!(spec.satisfies_concrete_theorem_gates());
        let eps = spec.epsilon().unwrap();
        assert!((eps - 1.0 / 18.0).abs() < 1e-12);
        // |W| = 1 gives a large kappa, above the 1/2 paper line.
        assert!(spec.kappa().unwrap() > 0.5);

        let dense = ThinProductSpec {
            n_outer: 64,
            middle_dim: 16,
            wanted: (0..64).flat_map(|i| (0..64).map(move |j| (i, j))).collect(),
        };
        assert!(!dense.in_concrete_regime());
        assert!(!dense.wanted_within_paper_bound());
        assert!(!dense.satisfies_concrete_theorem_gates());
    }

    #[test]
    fn cost_model_hooks_are_monotone_and_ordered() {
        // Larger D lowers the paper's offline bound at fixed N.
        let small_d = paper_ops_bound(1 << 20, 16).unwrap();
        let large_d = paper_ops_bound(1 << 20, 256).unwrap();
        assert!(large_d < small_d);
        // At the paper's max wanted size |W| = N^2/sqrt(D), the naive
        // baseline costs N^2 sqrt(D) while the paper quotes N^2/D^0.063;
        // the quoted bound wins by D^0.563 for every D > 1. Compare the
        // closed forms directly: materialising such a W is infeasible.
        let n: f64 = (1 << 20) as f64;
        for d in [4.0_f64, 16.0, 256.0, 4096.0] {
            let naive_max = n * n * d.sqrt();
            let bound = n * n / d.powf(0.063);
            assert!(
                bound < naive_max,
                "paper bound should beat max-size naive baseline at D={d}"
            );
        }
        // A tiny two-entry spec is below the regime where the asymptotic
        // bound bites; the naive cost there is just |W| * D exactly.
        let spec = ThinProductSpec {
            n_outer: 1 << 18,
            middle_dim: 4,
            wanted: vec![(0, 0), (1, 1)],
        };
        assert_eq!(naive_inner_product_cost(&spec), Some(8));
        // Query bound grows with D; graph bound grows with |W|.
        assert!(paper_query_bound(256).unwrap() > paper_query_bound(16).unwrap());
        let g1 = graph_time_bound(1000, 0.05, 100).unwrap();
        let g2 = graph_time_bound(1000, 0.05, 200).unwrap();
        assert!(g2 > g1);
    }

    #[test]
    fn triangle_counts_match_wanted_product_on_binary_inputs() {
        let spec = ThinProductSpec {
            n_outer: 3,
            middle_dim: 2,
            wanted: vec![(0, 0), (0, 2), (1, 1), (2, 0)],
        };
        // Left biadjacency (3x2), right biadjacency (2x3).
        let x01: Vec<u8> = vec![1, 0, 1, 1, 0, 1];
        let y01: Vec<u8> = vec![1, 0, 1, 0, 1, 1];
        let counts = count_triangles_through_pairs(&spec, &x01, &y01).unwrap();
        let x: Vec<i64> = x01.iter().map(|&b| b as i64).collect();
        let y: Vec<i64> = y01.iter().map(|&b| b as i64).collect();
        let prods = wanted_product_i64(&spec, &x, &y).unwrap();
        assert_eq!(counts, prods.iter().map(|&v| v as u64).collect::<Vec<_>>());
        // (0,0): k=0 only -> 1; (2,0): neither k -> 0.
        assert_eq!(counts[0], 1);
        assert_eq!(counts[3], 0);
        assert!(count_triangles_through_pairs(&spec, &[2, 0, 1, 1, 0, 1], &y01).is_err());
    }

    #[test]
    fn insertion_points_are_all_unimplemented_proposals() {
        for point in all_insertion_points() {
            let a = insertion_assessment(point);
            assert_eq!(a.point, point);
            assert!(!a.implemented);
            assert!(!a.blocking_difference.is_empty());
            assert!(!a.required_experiment.is_empty());
        }
    }
}
